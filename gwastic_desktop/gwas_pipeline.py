import subprocess
import pandas as pd
from gwastic_desktop.gwas_ai import GWASAI
import sys
import os

class GWAS:
    """GWAS class."""

    def __init__(self):
        self.gwas_ai = GWASAI()

    def vcf_to_bed(self, vcf_file, id_file, file_out, maf, geno):
        """Converts the vcf to bed files."""

        script_dir = os.path.dirname(__file__)  # <-- absolute dir the script is in



        if sys.platform.startswith('win'):
            #link_call = "windows/plink"
            rel_path = "windows/plink"
            abs_file_path = os.path.join(script_dir, rel_path)

        elif sys.platform.startswith('linux'):
            rel_path = "linux/plink"
            abs_file_path = os.path.join(script_dir, rel_path)


        print (abs_file_path)
        if id_file == None:
            process = subprocess.Popen([abs_file_path, "--vcf", vcf_file, "--make-bed", "--out", file_out,
                                        "--allow-extra-chr", "--set-missing-var-ids", "@:#", "--maf", maf,
                                        "--geno", geno, "--double-id"])

        else:
            process = subprocess.Popen([abs_file_path, "--vcf", vcf_file, "--make-bed", "--out", file_out,
                                        "--allow-extra-chr", "--set-missing-var-ids", "@:#", "--maf", maf,
                                        "--geno", geno, "--double-id", "--keep", id_file])


        process.wait()
        # read the log file
        #plink_log = open(file_out + '.log').read()
        plink_log = "Conversation done!"
        return plink_log

    def filter_out_missing(self, bed):
        import numpy as np
        import logging
        from pysnptools.util import log_in_place
        sid_batch_size = 1000  # number of SNPs to read at a time, decrease in case we run oput of memory
        # list of all NaN columns
        all_nan = []
        with log_in_place("read snp #", logging.INFO) as updater:
            for sid_start in range(0, bed.sid_count, sid_batch_size):
                updater(f"sid {sid_start:,} of {bed.sid_count:,}")
                # read a batch of SNPs
                snp_data_batch = bed[:, sid_start:sid_start + sid_batch_size].read()
                # find the columns that are all NaN (missing)
                nan_in_batch = np.isnan(snp_data_batch.val).all(axis=0)
                all_nan.append(nan_in_batch)

        # concatenate the list of all NaN columns into a single boolean array
        all_nan = np.concatenate(all_nan, axis=0)
        # print the number of all NaN columns
        logging.info(f"number of all missing columns is {np.sum(all_nan):,}. They are: ")
        logging.info(bed.sid[all_nan])
        # remove these SNPs from consideration
        bed_fixed = bed[:, ~all_nan]
        return bed_fixed

    def validate_gwas_input_files(self, bed_file, pheno_file):
        """Validate the input file.
        Phenotypic file must be ID_space_ID_value.
        fam file ids must match phenotype ids."""

        fam_data = open(bed_file.replace('.bed', '.fam'), 'r').readlines()
        fam_ids = []
        for line in fam_data:
            line = line.strip()
            line = line.split(' ')
            if len(line) > 1:
                fam_ids.append(line[0])
            else:
                return (False, "FAM file is not space delimited.")

        pheno_data = open(pheno_file, 'r').readlines()
        pheno_ids = []
        columns = None
        for line2 in pheno_data:
            line2 = line2.strip()
            line2 = line2.split(' ')
            columns = len(line2)
            if len(line2) > 1:
                pheno_ids.append(line2[0])
                #return (True, "Input files validated.")
            else:
                return (False, "Phenotpic file is not space delimited.")

        # Check whether all pheno ids are in the fam file
        check_ids = any(x in pheno_ids for x in fam_ids)

        if check_ids and columns == 3:
            return (True, "Input files validated.")
        elif columns != 3:
            return (False, "Invalid phenotypic file. File should contain 3 coloumns (FID IID Value)")
        else:
            return (False, "Phenotpic ID's does not match with .fam file IDs.")

    def start_gwas(self, bed_file, pheno_file ,chrom_mapping, algorithm, add_log, test_size, estimators,
                   gwas_result_name):
        from fastlmm.association import single_snp, single_snp_linreg
        from pysnptools.snpreader import Bed, Pheno
        import pysnptools.util as pstutil
        import time
        # First we have to validate the input files
        check_input_data = self.validate_gwas_input_files(bed_file, pheno_file)

        if check_input_data[0]:
            bed = Bed(str(bed_file), count_A1=False, chrom_map=chrom_mapping)
            pheno = Pheno(str(pheno_file))
            # replace original bed with one that has the same iids as the pheno
            bed, pheno = pstutil.intersect_apply([bed, pheno])
            bed_fixed = self.filter_out_missing(bed)

            # format numbers with commas and no decimals
            s3 = "Dataset after intersection:" + ' SNPs: ' + str(bed.sid_count) + ' Pheno IDs: ' + str(pheno.iid_count)
            add_log(s3, warn=True)
            # run single_snp with the fixed file
            add_log('Starting GWAS Analysis, this might take a while...')

            #self.plot_pheno(pheno.read().val)

            t1 = time.process_time()
            if algorithm == 'FaST-LMM':
                df = single_snp(bed_fixed, pheno, output_file_name=gwas_result_name)
                exchanged_dict = {v: k for k, v in chrom_mapping.items()}
                df['Chr'] = df['Chr'].replace(exchanged_dict)
            elif algorithm == 'Linear regression':
                df = single_snp_linreg(test_snps=bed_fixed, pheno=pheno, output_file_name=gwas_result_name)
                exchanged_dict = {v: k for k, v in chrom_mapping.items()}
                df['Chr'] = df['Chr'].replace(exchanged_dict)
            elif algorithm == 'Random Forest (AI)':
                df = pd.read_csv(bed_file.replace('bed', 'bim'), delimiter='\t')
                snp_ids = df.iloc[:, 1].tolist()
                df = self.gwas_ai.run_random_forest(bed_fixed.read().val, pheno.read().val, snp_ids, test_size,
                                                    estimators, gwas_result_name)
            elif algorithm == 'XGBoost (AI)':
                df = pd.read_csv(bed_file.replace('bed', 'bim'), delimiter='\t')
                snp_ids = df.iloc[:, 1].tolist()
                df = self.gwas_ai.run_xgboost(bed_fixed.read().val, pheno.read().val, snp_ids, test_size,
                                                    estimators, gwas_result_name)

            t2 = time.process_time()
            t3 = round((t2-t1)/ 60, 2)
            add_log('Final run time (minutes): ' + str(t3), 2)
            return df
        else:
            add_log(check_input_data[1], error=True)

    def plot_gwas(self, df, limit, algorithm, manhatten_plot_name, qq_plot_name):
        """Manhatten and qq-plot."""
        import matplotlib.pyplot as plt
        import geneview as gv

        if algorithm == 'FaST-LMM' or algorithm == 'Linear regression':
            df = df.sort_values(by=['Chr', 'ChrPos'])
            df['Chr'] = df['Chr'].astype(int)
            df['ChrPos'] = df['ChrPos'].astype(int)

            # common parameters for plotting
            plt_params = {
                "font.sans-serif": "Arial",
                "legend.fontsize": 14,
                "axes.titlesize": 18,
                "axes.labelsize": 16,
                "xtick.labelsize": 14,
                "ytick.labelsize": 14
            }
            plt.rcParams.update(plt_params)

            # Create a manhattan plot
            f, ax = plt.subplots(figsize=(12, 5), facecolor="w", edgecolor="k")
            #xtick = set(["chr" + i for i in list(map(str, chr_names))])
            _ = gv.manhattanplot(data=df,chrom='Chr', pos="ChrPos", pv="PValue", snp="SNP", marker=".", color=['#3B5488', '#d5536f'],
                              sign_marker_color="r", title="GWAS Manhatten Plot " + algorithm + '\n', #xtick_label_set=xtick,
                              xlabel="Chromosome", ylabel=r"$-log_{10}{(P)}$", sign_line_cols=["#D62728", "#2CA02C"],
                              hline_kws={"linestyle": "--", "lw": 1.3},   # 50000 bp
                                # hline_kws={"linestyle": "--", "lw": 1.3}, is_annotate_topsnp=True, ld_block_size=50000,sign_marker_p=1e-6,
                                 # 50000 bp

                                 text_kws={"fontsize": 12, "arrowprops": dict(arrowstyle="-", color="k", alpha=0.6)}, ax=ax)
            plt.tight_layout(pad=1)
            plt.savefig(manhatten_plot_name, dpi=100)

            # Create QQ plot
            f, ax = plt.subplots(figsize=(6, 6), facecolor="w", edgecolor="k")
            _ = gv.qqplot(data=df["PValue"], marker="o", title="GWAS QQ Plot " + algorithm + '\n',
                          xlabel=r"Expected $-log_{10}{(P)}$", ylabel=r"Observed $-log_{10}{(P)}$", ax=ax)
            #ax.set(ylim=(0, 20), xlim=(0, 20))
            plt.tight_layout(pad=1)
            plt.savefig(qq_plot_name, dpi=100)

        else:

            df['Chr'] = df['Chr'].astype(int)
            # chr_names = df['Chr'].unique()
            df['ChrPos'] = df['ChrPos'].astype(int)

            f, ax = plt.subplots(figsize=(12, 5), facecolor="w", edgecolor="k")

            _ = gv.manhattanplot(data=df, chrom='Chr', pos="ChrPos", pv="PValue", snp="SNP", logp=False,
                                  title="GWAS Manhatten Plot " + algorithm + '\n',color=['#3B5488', '#d5536f'],
                                  xlabel="Chromosome", ylabel=r"Feature Importance")

            plt.tight_layout(pad=1)
            plt.savefig(manhatten_plot_name, dpi=200)

    def plot_pheno(self, data):
        import matplotlib.pyplot as plt
        import seaborn as sns

        plt.hist(data, edgecolor='black')  # 'bins' defines the number of bins
        plt.xlabel('Phenotype Value')
        plt.ylabel('Frequency')
        plt.title('Phenotype Distribution')
        plt.show()









