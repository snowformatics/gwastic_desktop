#!/usr/bin/env python3
"""
gwastic_batch.py  –  Batch GWAS command-line tool powered by GWAStic Desktop
=============================================================================
Runs multiple GWAS analyses (one per phenotype .txt file) against a shared
set of PLINK BED/BIM/FAM files, using every selected algorithm.

Usage
-----
  python gwastic_batch.py --folder /path/to/data [OPTIONS]

Folder layout expected
----------------------
  <folder>/
    *.bed  (+ matching .bim and .fam)
    *.txt  (phenotype files, 3-column space-delimited: FID IID value)

Output
------
  <folder>/gwastic_results/<timestamp>_<phenotype>_<algorithm>/
    gwas_results.csv
    gwas_results_top10000.csv
    manhatten_plot_high.png
    qq_plot_high.png          (LMM / Linear regression only)
    log.txt

Dependencies
------------
  pip install fastlmm pysnptools xgboost scikit-learn geneview matplotlib pandas numpy
  (The gwastic_desktop package must also be importable – add its parent directory
   to PYTHONPATH or install it, e.g. pip install -e .)

Algorithms available
--------------------
  FaST-LMM            – Mixed linear model (best for population-structure correction)
  Linear regression   – Fast OLS regression
  Ridge Regression    – Penalised regression (alpha tunable via --ridge-alpha)
  Random Forest (AI)  – Ensemble feature-importance GWAS
  XGBoost (AI)        – Gradient-boosted feature-importance GWAS
"""

import argparse
import glob
import logging
import os
import shutil
import sys
import time
from datetime import datetime
from pathlib import Path

# ---------------------------------------------------------------------------
# Logging setup
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("gwastic_batch")

ALGORITHMS_ALL = [
    "FaST-LMM",
    "Linear regression",
    "Ridge Regression",
    "Random Forest (AI)",
    "XGBoost (AI)",
]


# ---------------------------------------------------------------------------
# Helper: simple add_log shim (the pipeline expects a callable)
# ---------------------------------------------------------------------------
def make_add_log(pheno_name: str, algorithm: str):
    prefix = f"[{pheno_name}|{algorithm}]"

    def _log(msg, error=False):
        if error:
            log.error(f"{prefix}  {msg}")
        else:
            log.info(f"{prefix}  {msg}")

    return _log


# ---------------------------------------------------------------------------
# Core batch runner
# ---------------------------------------------------------------------------
def run_batch(
    folder: str,
    algorithms: list,
    output_root: str,
    cov_file,  # str or None
    gb_goal: int,
    train_size: int,       # percentage, e.g. 70
    n_estimators: int,
    n_models: int,
    max_depth: int,
    n_jobs: int,
    aggregation: str,      # sum | mean | median
    ridge_alpha: float,
    snp_limit: str,
    skip_plots: bool,
    dry_run: bool,
):
    folder = os.path.abspath(folder)

    # ------------------------------------------------------------------
    # Locate BED file (expect exactly one per folder)
    # ------------------------------------------------------------------
    bed_files = glob.glob(os.path.join(folder, "*.bed"))
    if len(bed_files) == 0:
        log.error("No .bed file found in %s", folder)
        sys.exit(1)
    if len(bed_files) > 1:
        log.warning("Multiple .bed files found – using the first one: %s", bed_files[0])
    bed_path = bed_files[0]

    # ------------------------------------------------------------------
    # Locate phenotype files
    # ------------------------------------------------------------------
    pheno_files = glob.glob(os.path.join(folder, "*.txt"))
    if len(pheno_files) == 0:
        log.error("No phenotype .txt files found in %s", folder)
        sys.exit(1)

    log.info("BED file   : %s", bed_path)
    log.info("Phenotypes : %d file(s) found", len(pheno_files))
    log.info("Algorithms : %s", ", ".join(algorithms))

    if dry_run:
        log.info("=== DRY RUN – nothing will be executed ===")
        for pf in pheno_files:
            for alg in algorithms:
                log.info("  Would run  %-30s  with  %s", os.path.basename(pf), alg)
        return

    # ------------------------------------------------------------------
    # Lazy imports (heavy ML stack – only pay cost when actually running)
    # ------------------------------------------------------------------
    try:
        import numpy as np
        import matplotlib
        matplotlib.use("Agg")
        from pysnptools.snpreader import Bed, Pheno
        from gwastic_desktop.gwas_pipeline import GWAS
        from gwastic_desktop.helpers import HELPERS
    except ImportError as exc:
        log.error("Required package not found: %s", exc)
        log.error("Run:  pip install fastlmm pysnptools xgboost scikit-learn geneview matplotlib")
        log.error("And ensure the gwastic_desktop package is on PYTHONPATH.")
        sys.exit(1)

    gwas_obj = GWAS()
    helper   = HELPERS()

    # Build chromosome mapping from BIM once (shared across all runs)
    bim_path = bed_path.replace(".bed", ".bim")
    chrom_mapping = helper.replace_with_integers(bim_path)

    # ------------------------------------------------------------------
    # Validate every phenotype file before starting anything
    # ------------------------------------------------------------------
    log.info("Validating input files …")
    invalid = []
    for pf in pheno_files:
        ok, msg = gwas_obj.validate_gwas_input_files(bed_path, pf)
        if not ok:
            log.warning("  SKIP  %s  –  %s", os.path.basename(pf), msg)
            invalid.append(pf)
        else:
            log.info("  OK    %s", os.path.basename(pf))
    pheno_files = [p for p in pheno_files if p not in invalid]

    if len(pheno_files) == 0:
        log.error("No valid phenotype files remain after validation. Exiting.")
        sys.exit(1)

    # ------------------------------------------------------------------
    # Load & fix BED (once – reuse across all phenotype × algorithm runs)
    # ------------------------------------------------------------------
    log.info("Loading BED file and filtering fully-missing SNPs …")
    bed_raw   = Bed(str(bed_path), count_A1=False, chrom_map=chrom_mapping)
    bed_fixed = gwas_obj.filter_out_missing(bed_raw)
    log.info("  SNPs after filtering: %d", bed_fixed.sid_count)

    # Covariate object (optional, shared)
    cov_obj = None
    if cov_file:
        if os.path.isfile(cov_file):
            cov_obj = Pheno(str(cov_file))
            log.info("Covariate file loaded: %s", cov_file)
        else:
            log.warning("Covariate file not found – running without covariates: %s", cov_file)

    # test-size fraction from training percentage
    test_size = (100 - train_size) / 100.0

    # ------------------------------------------------------------------
    # Main loop: phenotype × algorithm
    # ------------------------------------------------------------------
    total      = len(pheno_files) * len(algorithms)
    done       = 0
    t_global   = time.time()

    for pheno_path in sorted(pheno_files):
        pheno_name = Path(pheno_path).stem          # e.g. "Browning_WGS_SAMEA_GWAS"
        pheno_obj  = Pheno(str(pheno_path))

        for algorithm in algorithms:
            done += 1
            log.info("")
            log.info("=" * 70)
            log.info("Run %d / %d  |  Phenotype: %s  |  Algorithm: %s",
                     done, total, pheno_name, algorithm)
            log.info("=" * 70)

            add_log = make_add_log(pheno_name, algorithm)

            # ---- per-run output directory --------------------------------
            ts_tag = datetime.now().strftime("%Y%m%d_%H%M%S")
            alg_tag = algorithm.replace(" ", "_").replace("(", "").replace(")", "")
            run_dir = os.path.join(
                output_root,
                f"{ts_tag}_{pheno_name}_{alg_tag}",
            )
            os.makedirs(run_dir, exist_ok=True)

            # ---- temp file names (written to CWD then copied) -----------
            gwas_result_name     = "gwas_results.csv"
            gwas_result_name_top = "gwas_results_top10000.csv"
            manhattan_plot       = "manhatten_plot.png"
            manhattan_plot_high  = "manhatten_plot_high.png"
            qq_plot              = "qq_plot.png"
            qq_plot_high         = "qq_plot_high.png"

            t_run = time.time()
            try:
                # ---- Run the selected algorithm --------------------------
                if algorithm in ("FaST-LMM", "Linear regression"):
                    gwas_df, df_plot = gwas_obj.run_gwas_lmm(
                        bed_fixed, pheno_obj, chrom_mapping, add_log,
                        gwas_result_name, algorithm, bed_path,
                        cov_obj, gb_goal,
                    )

                elif algorithm == "Ridge Regression":
                    gwas_df, df_plot = gwas_obj.run_gwas_ridge(
                        bed_fixed, pheno_obj, bed_path, test_size,
                        ridge_alpha, gwas_result_name, chrom_mapping,
                        add_log, n_models, aggregation,
                    )

                elif algorithm == "Random Forest (AI)":
                    gwas_df, df_plot = gwas_obj.run_gwas_rf(
                        bed_fixed, pheno_obj, bed_path, test_size,
                        n_estimators, gwas_result_name, chrom_mapping,
                        add_log, n_models, n_jobs, aggregation,
                    )

                elif algorithm == "XGBoost (AI)":
                    gwas_df, df_plot = gwas_obj.run_gwas_xg(
                        bed_fixed, pheno_obj, bed_path, test_size,
                        n_estimators, gwas_result_name, chrom_mapping,
                        add_log, n_models, max_depth, n_jobs, aggregation,
                    )

                else:
                    log.error("Unknown algorithm '%s' – skipping.", algorithm)
                    continue

                # ---- Plots ----------------------------------------------
                if not skip_plots:
                    try:
                        gwas_obj.plot_gwas(
                            df_plot, snp_limit, algorithm,
                            manhattan_plot, qq_plot, chrom_mapping,
                        )
                        add_log("Plots generated.")
                    except Exception as e:
                        add_log(f"Plot generation failed (non-fatal): {e}", error=True)

                # ---- Top-10000 file -------------------------------------
                try:
                    import pandas as pd
                    df_top = pd.read_csv(gwas_result_name).head(10000)
                    df_top.to_csv(gwas_result_name_top, index=False)
                except Exception:
                    pass

                # ---- Copy all output files to run_dir -------------------
                elapsed = round((time.time() - t_run) / 60, 2)
                _copy_results(
                    run_dir,
                    gwas_result_name, gwas_result_name_top,
                    manhattan_plot_high, qq_plot_high,
                    algorithm, pheno_name, bed_path, pheno_path,
                    test_size, n_estimators, n_models, max_depth,
                    aggregation, elapsed,
                )
                log.info("Results saved to: %s", run_dir)

            except Exception as exc:
                log.error("FAILED  %s  with  %s:  %s", pheno_name, algorithm, exc, exc_info=True)
                _write_error_log(run_dir, pheno_name, algorithm, str(exc))

            finally:
                # Clean up temp files in CWD
                for tmp in [gwas_result_name, gwas_result_name_top,
                            manhattan_plot, manhattan_plot_high,
                            qq_plot, qq_plot_high]:
                    if os.path.exists(tmp):
                        os.remove(tmp)

    elapsed_total = round((time.time() - t_global) / 60, 2)
    log.info("")
    log.info("All runs complete. Total time: %s min", elapsed_total)
    log.info("Results directory: %s", output_root)


# ---------------------------------------------------------------------------
# Helper: copy results + write log
# ---------------------------------------------------------------------------
def _copy_results(
    run_dir, gwas_result_name, gwas_result_name_top,
    manhattan_plot_high, qq_plot_high,
    algorithm, pheno_name, bed_path, pheno_path,
    test_size, n_estimators, n_models, max_depth,
    aggregation, elapsed_min,
):
    for src in [gwas_result_name, gwas_result_name_top,
                manhattan_plot_high, qq_plot_high]:
        if os.path.exists(src):
            shutil.copy(src, os.path.join(run_dir, src))

    with open(os.path.join(run_dir, "log.txt"), "w") as f:
        f.write(f"Algorithm        : {algorithm}\n")
        f.write(f"BED file         : {bed_path}\n")
        f.write(f"Phenotype file   : {pheno_path}\n")
        f.write(f"Phenotype name   : {pheno_name}\n")
        f.write(f"Test size        : {round(test_size * 100)}%\n")
        f.write(f"Nr. of trees     : {n_estimators}\n")
        f.write(f"Nr. of models    : {n_models}\n")
        f.write(f"Max depth        : {max_depth}\n")
        f.write(f"Aggregation      : {aggregation}\n")
        f.write(f"Runtime (min)    : {elapsed_min}\n")


def _write_error_log(run_dir, pheno_name, algorithm, error_msg):
    with open(os.path.join(run_dir, "ERROR.txt"), "w") as f:
        f.write(f"Phenotype : {pheno_name}\n")
        f.write(f"Algorithm : {algorithm}\n")
        f.write(f"Error     : {error_msg}\n")


# ---------------------------------------------------------------------------
# CLI argument parsing
# ---------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "GWAStic Batch – run multiple GWAS analyses from a folder.\n\n"
            "All *.txt files in --folder are treated as phenotype files.\n"
            "The *.bed file (+ .bim / .fam) must also be in the same folder."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # --- Required ---------------------------------------------------------
    p.add_argument(
        "--folder", "-f", required=True,
        help="Folder containing the .bed file and phenotype .txt files.",
    )

    # --- Algorithm selection ----------------------------------------------
    alg_group = p.add_mutually_exclusive_group()
    alg_group.add_argument(
        "--algorithms", "-a", nargs="+",
        choices=ALGORITHMS_ALL,
        metavar="ALG",
        help=(
            "One or more algorithms to run (space-separated).\n"
            "Choices: " + " | ".join(ALGORITHMS_ALL) + "\n"
            "Default: ALL algorithms are run."
        ),
    )
    alg_group.add_argument(
        "--all-algorithms", action="store_true", default=False,
        help="Run all 5 algorithms (default behaviour).",
    )

    # --- Output -----------------------------------------------------------
    p.add_argument(
        "--output", "-o",
        default=None,
        help=(
            "Root directory for results. "
            "Defaults to <folder>/gwastic_results/."
        ),
    )

    # --- Optional inputs --------------------------------------------------
    p.add_argument(
        "--cov", default=None,
        help="Path to an optional covariate file (3-column: FID IID value).",
    )

    # --- LMM parameters ---------------------------------------------------
    lmm = p.add_argument_group("FaST-LMM / Linear regression options")
    lmm.add_argument(
        "--gb-goal", type=int, default=0,
        help="Memory goal in GB for FaST-LMM (0 = auto, default: 0).",
    )
    lmm.add_argument(
        "--snp-limit", default="",
        help="Limit plot to top N SNPs (leave empty for all, default: all).",
    )

    # --- ML parameters ----------------------------------------------------
    ml = p.add_argument_group("ML algorithm options (RF, XGBoost, Ridge)")
    ml.add_argument(
        "--train-size", type=int, default=70, metavar="PCT",
        help="Training set size in percent (default: 70, i.e. 30%% test).",
    )
    ml.add_argument(
        "--n-estimators", type=int, default=200,
        help="Number of trees for RF / XGBoost (default: 200).",
    )
    ml.add_argument(
        "--n-models", type=int, default=1,
        help="Number of model iterations; SNP effects are aggregated (default: 1).",
    )
    ml.add_argument(
        "--max-depth", type=int, default=3,
        help="Max tree depth for XGBoost (default: 3).",
    )
    ml.add_argument(
        "--n-jobs", type=int, default=-1,
        help="Parallel jobs for RF / XGBoost (-1 = all CPUs, default: -1).",
    )
    ml.add_argument(
        "--aggregation", choices=["sum", "mean", "median"], default="sum",
        help="How to aggregate SNP effects across model iterations (default: sum).",
    )
    ml.add_argument(
        "--ridge-alpha", type=float, default=1.0,
        help="Regularisation strength for Ridge Regression (default: 1.0).",
    )

    # --- Misc -------------------------------------------------------------
    p.add_argument(
        "--skip-plots", action="store_true",
        help="Skip Manhattan and QQ plot generation (faster).",
    )
    p.add_argument(
        "--dry-run", action="store_true",
        help="Print what would be run without actually running anything.",
    )
    p.add_argument(
        "--list-algorithms", action="store_true",
        help="Print available algorithms and exit.",
    )

    return p


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
def main():
    parser = build_parser()
    args   = parser.parse_args()

    if args.list_algorithms:
        print("\nAvailable GWAS algorithms:")
        for i, alg in enumerate(ALGORITHMS_ALL, 1):
            print(f"  {i}. {alg}")
        print()
        sys.exit(0)

    # Resolve algorithm list
    if args.algorithms:
        algorithms = args.algorithms
    else:
        algorithms = ALGORITHMS_ALL  # default: run all

    # Resolve output directory
    output_root = args.output or os.path.join(args.folder, "gwastic_results")
    os.makedirs(output_root, exist_ok=True)

    # Validate train-size
    if not (1 <= args.train_size <= 99):
        parser.error("--train-size must be between 1 and 99.")

    log.info("GWAStic Batch CLI")
    log.info("  Folder          : %s", args.folder)
    log.info("  Output root     : %s", output_root)
    log.info("  Algorithms      : %s", ", ".join(algorithms))
    log.info("  Covariate file  : %s", args.cov or "none")
    log.info("  GB goal         : %s", args.gb_goal)
    log.info("  Train size      : %d%%", args.train_size)
    log.info("  Estimators      : %d", args.n_estimators)
    log.info("  Models          : %d", args.n_models)
    log.info("  Max depth       : %d", args.max_depth)
    log.info("  Jobs            : %d", args.n_jobs)
    log.info("  Aggregation     : %s", args.aggregation)
    log.info("  Ridge alpha     : %.4f", args.ridge_alpha)
    log.info("  SNP limit       : %s", args.snp_limit or "all")
    log.info("  Skip plots      : %s", args.skip_plots)
    log.info("  Dry run         : %s", args.dry_run)

    run_batch(
        folder        = args.folder,
        algorithms    = algorithms,
        output_root   = output_root,
        cov_file      = args.cov,
        gb_goal       = args.gb_goal,
        train_size    = args.train_size,
        n_estimators  = args.n_estimators,
        n_models      = args.n_models,
        max_depth     = args.max_depth,
        n_jobs        = args.n_jobs,
        aggregation   = args.aggregation,
        ridge_alpha   = args.ridge_alpha,
        snp_limit     = args.snp_limit,
        skip_plots    = args.skip_plots,
        dry_run       = args.dry_run,
    )


if __name__ == "__main__":
    main()