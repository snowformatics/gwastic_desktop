import subprocess
import pandas as pd
from gwastic_desktop.gwas_ai import GWASAI
import sys
import os
from gwastic_desktop.helpers import HELPERS


class GWAS:
    """GWAS class."""

    def __init__(self):
        self.gwas_ai = GWASAI()
        self.helper = HELPERS()

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
            #os.chmod(abs_file_path, 0o755)

        #print (abs_file_path)
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
            return True, "Input files validated."
        elif columns != 3:
            return (False, "Invalid phenotypic file. File should contain 3 coloumns (FID IID Value)")
        else:
            return (False, "Phenotpic ID's does not match with .fam file IDs.")

    def start_gwas(self, bed_file, pheno_file, chrom_mapping, algorithm, add_log, test_size, model_nr, estimators, leave_chr_set,
                   max_dep_set, gwas_result_name, genomic_predict, genomic_predict_name):
        from fastlmm.association import single_snp, single_snp_linreg
        from pysnptools.snpreader import Bed, Pheno
        import pysnptools.util as pstutil
        import time

        print (test_size, estimators, model_nr)
        # First we have to validate the input files
        check_input_data = self.validate_gwas_input_files(bed_file, pheno_file)

        if check_input_data[0]:
            t1 = time.time()

            bed = Bed(str(bed_file), count_A1=False, chrom_map=chrom_mapping)
            pheno = Pheno(str(pheno_file))

            if genomic_predict:
                bed_gp = Bed(str(bed_file), count_A1=False, chrom_map=chrom_mapping)
                pheno_gp = Pheno(str(pheno_file))
            else:
                bed_gp = None
                pheno_gp = None

            # replace original bed with one that has the same iids as the pheno
            bed, pheno = pstutil.intersect_apply([bed, pheno])
            bed_fixed = self.filter_out_missing(bed)

            # format numbers with commas and no decimals
            s3 = "Dataset after intersection:" + ' SNPs: ' + str(bed.sid_count) + ' Pheno IDs: ' + str(pheno.iid_count)
            add_log(s3, warn=True)
            # run single_snp with the fixed file
            add_log('Starting Analysis, this might take a while...')
            #print (algorithm)

            if algorithm == 'FaST-LMM':
                df = single_snp(bed_fixed, pheno, leave_out_one_chrom=leave_chr_set, output_file_name=gwas_result_name)
                exchanged_dict = {v: k for k, v in chrom_mapping.items()}
                df['Chr'] = df['Chr'].replace(exchanged_dict)
            elif algorithm == 'Linear regression':
                df = single_snp_linreg(test_snps=bed_fixed, pheno=pheno, output_file_name=gwas_result_name)
                exchanged_dict = {v: k for k, v in chrom_mapping.items()}
                df['Chr'] = df['Chr'].replace(exchanged_dict)

            elif algorithm == 'Random Forest (AI)':
                dataframes = []
                for i in range(int(model_nr)):
                    #print(i)
                    df = pd.read_csv(bed_file.replace('bed', 'bim'), delimiter='\t')
                    snp_ids = df.iloc[:, 1].tolist()
                    df = self.gwas_ai.run_random_forest(bed_fixed.read().val, pheno.read().val, snp_ids, test_size,
                                                  estimators, gwas_result_name, bed_gp, pheno_gp, genomic_predict,
                                                  genomic_predict_name, model_nr)
                    dataframes.append(df)
                # We merge the models and calculate the sum effect
                if genomic_predict:
                    df = df
                else:
                    df = self.helper.merge_models(dataframes)
            elif algorithm == 'XGBoost (AI)':
                dataframes = []
                for i in range(model_nr):
                    #print (i)
                    df = pd.read_csv(bed_file.replace('bed', 'bim'), delimiter='\t')
                    snp_ids = df.iloc[:, 1].tolist()
                    df = self.gwas_ai.run_xgboost(bed_fixed.read().val, pheno.read().val, snp_ids, test_size,
                                                  estimators, str(i) + gwas_result_name, bed_gp, pheno_gp, genomic_predict,
                                                  genomic_predict_name, max_dep_set, model_nr)
                    dataframes.append(df)
                if genomic_predict:
                    df = df
                else:
                    df = self.helper.merge_models(dataframes)
            elif algorithm == 'GP_LMM':
                #dataframes = []
                print (algorithm)
                for i in range(model_nr):
                    #print (i)
                    df = pd.read_csv(bed_file.replace('bed', 'bim'), delimiter='\t')
                    snp_ids = df.iloc[:, 1].tolist()
                    #df = self.gwas_ai.run_xgboost(bed_fixed.read().val, pheno.read().val, snp_ids, test_size,
                                                 # estimators, str(i) + gwas_result_name, bed_gp, pheno_gp, genomic_predict,
                                                 # genomic_predict_name, max_dep_set, model_nr)
                    df = self.gwas_ai.run_lmm_gp(bed_fixed.read().val, pheno.read().val, snp_ids, test_size,
                                                  estimators, str(i) + gwas_result_name, bed_gp, pheno_gp, genomic_predict,
                                                  genomic_predict_name, max_dep_set, model_nr)
                    #dataframes.append(df)
                if genomic_predict:
                    df = df
                #else:
                    #df = self.helper.merge_models(dataframes)

            t2 = time.time()
            t3 = round((t2-t1)/ 60, 2)
            df.to_csv(gwas_result_name, index=0)
            add_log('Final run time (minutes): ' + str(t3))
            return df
        else:
            add_log(check_input_data[1], error=True)

    def plot_gwas(self, df, limit, algorithm, manhatten_plot_name, qq_plot_name):
        """Manhatten and qq-plot."""
        import matplotlib.pyplot as plt
        import geneview as gv

        # Extract the top10 SNPs and use the value as significant marker label threshold

        if algorithm == 'FaST-LMM' or algorithm == 'Linear regression':
            # sign_limit = df.nsmallest(125, 'PValue')
            # sign_limit = sign_limit.tail(1)
            # sign_limit = float(sign_limit['PValue'])
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
            _ = gv.manhattanplot(data=df,chrom='Chr', pos="ChrPos", pv="PValue", snp="SNP", marker=".",color=['#4297d8', '#eec03c','#423496','#495227','#d50b6f','#e76519','#d580b7','#84d3ac'],
                              sign_marker_color="r", title="GWAS Manhatten Plot " + algorithm + '\n', #xtick_label_set=xtick,
                              xlabel="Chromosome", ylabel=r"$-log_{10}{(P)}$", sign_line_cols=["#D62728", "#2CA02C"],
                              hline_kws={"linestyle": "--", "lw": 1.3},#, sign_marker_p=1e-9, is_annotate_topsnp=True,
                                 text_kws={"fontsize": 12, "arrowprops": dict(arrowstyle="-", color="k", alpha=0.6)}, ax=ax)
            plt.tight_layout(pad=1)
            plt.savefig(manhatten_plot_name)

            # Create QQ plot
            f, ax = plt.subplots(figsize=(6, 6), facecolor="w", edgecolor="k")
            _ = gv.qqplot(data=df["PValue"], marker="o", title="GWAS QQ Plot " + algorithm + '\n',
                          xlabel=r"Expected $-log_{10}{(P)}$", ylabel=r"Observed $-log_{10}{(P)}$", ax=ax)
            #ax.set(ylim=(0, 20), xlim=(0, 20))
            plt.tight_layout(pad=1)
            plt.savefig(qq_plot_name)

        else:
            # sign_limit = df.nlargest(5, 'PValue')
            # print(sign_limit)
            # sign_limit = sign_limit.tail(1)
            # sign_limit = float(sign_limit['PValue'])
            # print (sign_limit)
            df['Chr'] = df['Chr'].astype(int)
            df['ChrPos'] = df['ChrPos'].astype(int)
            #print(df)
            #print (type(df))
            #print(df.dtypes)
            #print (df.info())

            f, ax = plt.subplots(figsize=(12, 5), facecolor="w", edgecolor="k")
            _ = gv.manhattanplot(data=df, chrom='Chr', pos="ChrPos", pv="PValue", snp="SNP", logp=False,
                                  title="GWAS Manhatten Plot " + algorithm + '\n',color=['#4297d8', '#eec03c','#423496','#495227','#d50b6f','#e76519','#d580b7','#84d3ac'],
                                  xlabel="Chromosome", ylabel=r"Feature Importance)")#, sign_marker_p=sign_limit, is_annotate_topsnp=True)
            plt.tight_layout(pad=1)
            plt.savefig(manhatten_plot_name)

    def plot_pheno(self, data):
        import matplotlib.pyplot as plt

        plt.hist(data, edgecolor='black')  # 'bins' defines the number of bins
        plt.xlabel('Phenotype Value')
        plt.ylabel('Frequency')
        plt.title('Phenotype Distribution')
        plt.show()

    def plot_gp(self, data, gp_plot_name, algorithm):
        """Bland-Altman Plot for the real and predicted phenotype values."""
        import matplotlib.pyplot as plt
        import numpy as np

        # Calculate means and differences
        means = (data['Predicted_Value'] + data['Pheno_Value']) / 2
        differences = data['Predicted_Value'] - data['Pheno_Value']
        mean_difference = np.mean(differences)
        std_difference = np.std(differences)

        # Plotting the Bland-Altman plot
        plt.figure(figsize=(10, 6))
        plt.scatter(means, differences, color='blue')
        plt.axhline(mean_difference, color='red', linestyle='--', label='Mean Difference')
        plt.axhline(mean_difference + 1.96 * std_difference, color='green', linestyle='--', label='Upper Limit of Agreement')
        plt.axhline(mean_difference - 1.96 * std_difference, color='green', linestyle='--', label='Lower Limit of Agreement')
        plt.xlabel('Mean Value')
        plt.ylabel('Difference')
        plt.title('Bland-Altman Plot (' + algorithm + ')')
        plt.legend()
        plt.tight_layout(pad=1)
        plt.savefig(gp_plot_name)
        #plt.show()

        # data = data.sort_values(by='Difference', ascending=False)
        # # Plotting
        # plt.figure(figsize=(10, 6))
        # # Scatter plot of predicted values
        # plt.scatter(data['ID1'], data['Predicted_Value'], label='Predicted Values')
        # # Adding error bars
        # #plt.errorbar(range(len(data['Predicted_Value'])), data['Predicted_Value'], yerr=data['Difference'], fmt='o',
        #              #ecolor='red', capsize=5)
        # # Optional: plot real values
        # plt.scatter(data['ID1'], data['Pheno_Value'], color='green', label='Real Values')
        # # Labels and Title
        # plt.xlabel('Sample IDs')
        # plt.ylabel('Values')
        # plt.title('Genomic Prediction Plot')
        # plt.xticks(rotation=90)
        # plt.legend()
        # #plt.grid(True)
        # plt.show()









