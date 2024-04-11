import subprocess
import sys
import os
from gwastic_desktop.helpers import HELPERS
from gwastic_desktop.gwas_ai import GWASAI
from fastlmm.association import single_snp, single_snp_linreg
from pysnptools.snpreader import Bed, Pheno
import pysnptools.util as pstutil
import time
from fastlmm.inference import FastLMM
from sklearn.model_selection import KFold
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import xgboost as xgb
from sklearn.ensemble import RandomForestRegressor
import geneview as gv
import numpy as np
import logging
from pysnptools.util import log_in_place
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import seaborn as sns

plt.switch_backend('Agg')

class GWAS:
    """GWAS class."""

    def __init__(self):
        self.gwas_ai = GWASAI()
        self.helper = HELPERS()

    def vcf_to_bed(self, vcf_file, id_file, file_out, maf, geno):
        """Converts the vcf to bed files."""

        script_dir = os.path.dirname(__file__)  # <-- absolute dir the script is in

        if sys.platform.startswith('win'):
            rel_path = "windows/plink"
            abs_file_path = os.path.join(script_dir, rel_path)
        elif sys.platform.startswith('linux'):
            rel_path = "linux/plink"
            abs_file_path = os.path.join(script_dir, rel_path)
            os.chmod(abs_file_path, 0o755)
        elif sys.platform.startswith('darwin'):
            rel_path = "mac/plink"
            abs_file_path = os.path.join(script_dir, rel_path)
            os.chmod(abs_file_path, 0o755)

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

    def validate_gwas_input_files(self, bed_file, pheno_file, ):
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


    def run_gwas_lmm(self, bed_fixed, pheno, chrom_mapping, add_log, gwas_result_name, algorithm):
        """GWAS using LMM and linear regression methods from fast-lmm library."""
        t1 = time.time()
        if algorithm == 'FaST-LMM':
            df_lmm_gwas = single_snp(bed_fixed, pheno, output_file_name=gwas_result_name)

        elif algorithm == 'Linear regression':
            df_lmm_gwas = single_snp_linreg(test_snps=bed_fixed, pheno=pheno, output_file_name=gwas_result_name)
        #print (df_lmm_gwas)
        df_lmm_gwas.dropna(subset=['PValue'], inplace=True)
        # we create one df for the plotting with ints as chr
        df_plot = df_lmm_gwas.copy(deep=True)

        if len(chrom_mapping) > 0:
            reversed_chrom_map = {value: key for key, value in chrom_mapping.items()}
            df_lmm_gwas["Chr"] = df_lmm_gwas["Chr"].apply(lambda x: reversed_chrom_map[x])

        df_lmm_gwas.to_csv(gwas_result_name, index=0)
        t2 = time.time()
        t3 = round((t2 - t1) / 60, 2)
        add_log('Final run time (minutes): ' + str(t3))
        return df_lmm_gwas, df_plot


    def run_gwas_xg(self, bed_fixed, pheno, bed_file, test_size, estimators, gwas_result_name, chrom_mapping, add_log,
                    model_nr, max_dep_set):
        """GWAS using Random Forest with cross validation."""

        t1 = time.time()
        dataframes = []
        df_bim = pd.read_csv(bed_file.replace('bed', 'bim'), delimiter='\t', header=None)
        df_bim.columns = ['Chr', 'SNP', 'NA', 'ChrPos', 'NA', 'NA']

        snp_data = bed_fixed.read().val
        snp_data[np.isnan(snp_data)] = -1

        for i in range(int(model_nr)):
            add_log('Model Iteration: ' + str(i + 1))
           # bed_fixed.read().val[np.isnan(bed_fixed.read().val)] = -1

            # Split data into training and testing sets
            X_train, X_test, y_train, y_test = train_test_split(snp_data, pheno.read().val,
                                                                test_size=test_size)

            scaler = StandardScaler()
            X_train = scaler.fit_transform(X_train)
            xg_model = xgb.XGBRegressor(n_estimators=estimators, learning_rate=0.1, max_depth=max_dep_set)
            xg_model.fit(X_train, y_train)

            data = []
            snp_ids = df_bim.iloc[:, 1].tolist()
            for col, score in zip(snp_ids, xg_model.feature_importances_):
                data.append((col, score))
            # Convert the list of tuples into a DataFrame
            df = pd.DataFrame(data, columns=['SNP', 'PValue'])
            df = pd.merge(df, df_bim, on='SNP', how='left')
            del df['NA']
            df['Chr'] = df['Chr'].replace(chrom_mapping)
            dataframes.append(df)

        df_all_xg = self.helper.merge_models(dataframes)
        df_all_xg = df_all_xg.sort_values(by='PValue', ascending=False)
        # we create one df for the plotting with ints as chr
        df_plot = df_all_xg.copy(deep=True)

        if len(chrom_mapping) > 0:
            reversed_chrom_map = {value: key for key, value in chrom_mapping.items()}
            df_all_xg["Chr"] = df_all_xg["Chr"].apply(lambda x: reversed_chrom_map[x])

        df_all_xg.to_csv(gwas_result_name, index=0)
        t2 = time.time()
        t3 = round((t2 - t1) / 60, 2)
        add_log('Final run time (minutes): ' + str(t3))
        return df_all_xg, df_plot

    def run_gwas_rf(self, bed_fixed, pheno, bed_file, test_size, estimators, gwas_result_name, chrom_mapping, add_log,
                    model_nr):
        """GWAS using Random Forest with cross validation."""
        t1 = time.time()
        dataframes = []
        df_bim = pd.read_csv(bed_file.replace('bed', 'bim'), delimiter='\t', header=None)
        df_bim.columns = ['Chr', 'SNP', 'NA', 'ChrPos', 'NA', 'NA']

        snp_data = bed_fixed.read().val
        snp_data[np.isnan(snp_data)] = -1

        for i in range(int(model_nr)):
            add_log('Model Iteration: ' + str(i + 1))
            # Split data into training and testing sets
            X_train, X_test, y_train, y_test = train_test_split(snp_data, pheno.read().val,
                                                                test_size=test_size)
            scaler = StandardScaler()
            X_train = scaler.fit_transform(X_train)
            rf_model = RandomForestRegressor(n_estimators=estimators, n_jobs=-1)
            rf_model.fit(X_train, y_train.ravel())
            data = []
            snp_ids = df_bim.iloc[:, 1].tolist()
            for col, score in zip(snp_ids, rf_model.feature_importances_):
                data.append((col, score))
            # Convert the list of tuples into a DataFrame
            df = pd.DataFrame(data, columns=['SNP', 'PValue'])
            df = pd.merge(df, df_bim, on='SNP', how='left')
            del df['NA']
            df['Chr'] = df['Chr'].replace(chrom_mapping)
            dataframes.append(df)

        df_all_rf = self.helper.merge_models(dataframes)
        df_all_rf = df_all_rf.sort_values(by='PValue', ascending=False)
        # we create one df for the plotting with ints as chr
        df_plot = df_all_rf.copy(deep=True)

        if len(chrom_mapping) > 0:
            reversed_chrom_map = {value: key for key, value in chrom_mapping.items()}
            df_all_rf["Chr"] = df_all_rf["Chr"].apply(lambda x: reversed_chrom_map[x])

        df_all_rf.to_csv(gwas_result_name, index=0)

        t2 = time.time()
        t3 = round((t2 - t1) / 60, 2)
        add_log('Final run time (minutes): ' + str(t3))
        return df_all_rf, df_plot

    def run_lmm_gp(self, bed_fixed, pheno, genomic_predict_name, model_nr, add_log, bed_file, chrom_mapping):
        """Genomic Prediction using LMM predictor from fast-lmm library."""

        bed_data = Bed(str(bed_file), count_A1=False, chrom_map=chrom_mapping)

        kf = KFold(n_splits=model_nr + 1)
        results = []
        i = 1
        # Start cross-validation
        for train_indices, test_indices in kf.split(range(bed_fixed.iid_count)):
            add_log('Model Iteration: ' + str(i))
            # Split data into training and testing
            train = bed_fixed[train_indices, :]
            #test = bed_fixed[test_indices, :]
            # Train the model
            fastlmm = FastLMM(GB_goal=2)
            fastlmm.fit(K0_train=train, y=pheno)
            # Test the model
            mean, covariance = fastlmm.predict(K0_whole_test=bed_data)
            df = pd.DataFrame({'Column1': mean.val[:, 0], 'Column2': mean.iid[:, 0]})

            #bed_data = bed_fixed.iid
            bed_data2 = bed_data.iid
            pheno_data2 = pheno.iid
            df_gp = pd.DataFrame(bed_data2, columns=['ID1', 'BED_ID2'])
            predicted_values = df['Column1']
            df_gp['Predicted_Value'] = predicted_values

            df_pheno = pd.DataFrame(pheno_data2, columns=['ID1', 'Pheno_ID2'])
            df_pheno['Pheno_Value'] = pheno.read().val
            merged_df = df_gp.merge(df_pheno, on='ID1', how='outer')
            merged_df['Predicted_Value'] = merged_df['Predicted_Value'].astype(float)
            merged_df['Predicted_Value'] = merged_df['Predicted_Value'].round(3)

            results.append(merged_df)
            i += 1

        df_all = self.helper.merge_gp_models(results)
        df_all.to_csv(genomic_predict_name, index=False)
        return df_all


    def run_gp_rf(self, bed_fixed, pheno, bed_file, test_size, estimators, genomic_predict_name, chrom_mapping, add_log,
                    model_nr):
        """Genomic Prediction using Random Forest with cross validation."""
        t1 = time.time()
        dataframes = []
        df_bim = pd.read_csv(bed_file.replace('bed', 'bim'), delimiter='\t', header=None)
        df_bim.columns = ['Chr', 'SNP', 'NA', 'ChrPos', 'NA', 'NA']

        # training data
        snp_data = bed_fixed.read().val
        snp_data[np.isnan(snp_data)] = -1
        # entire data for training
        bed_data = Bed(str(bed_file), count_A1=False, chrom_map=chrom_mapping)
        snp_data_all = bed_data.read().val
        snp_data_all[np.isnan(snp_data_all)] = -1
        bed_data2 = bed_data.iid
        pheno_data2 = pheno.iid

        for i in range(int(model_nr)):
            add_log('Model Iteration: ' + str(i + 1))
            # Split data into training and testing sets
            X_train, X_test, y_train, y_test = train_test_split(snp_data, pheno.read().val,
                                                                test_size=test_size)

            #scaler = StandardScaler()
            #X_train = scaler.fit_transform(X_train)
            rf_model = RandomForestRegressor(n_estimators=estimators)

            rf_model.fit(X_train, y_train.ravel())
            predicted_values = rf_model.predict(snp_data_all)

            df_gp = pd.DataFrame(bed_data2, columns=['ID1', 'BED_ID2'])
            df_gp['Predicted_Value'] = predicted_values

            df_pheno = pd.DataFrame(pheno_data2, columns=['ID1', 'Pheno_ID2'])
            df_pheno['Pheno_Value'] = pheno.read().val
            merged_df = df_gp.merge(df_pheno, on='ID1', how='outer')
            merged_df['Predicted_Value'] = merged_df['Predicted_Value'].astype(float)
            merged_df['Predicted_Value'] = merged_df['Predicted_Value'].round(5)
            dataframes.append(merged_df)


        df_all = self.helper.merge_gp_models(dataframes)
        df_all.to_csv(genomic_predict_name, index=False)
        t2 = time.time()
        t3 = round((t2 - t1) / 60, 2)
        add_log('Final run time (minutes): ' + str(t3))
        return df_all

    def run_gp_xg(self, bed_fixed, pheno, bed_file, test_size, estimators, genomic_predict_name, chrom_mapping, add_log,
                    model_nr, max_dep_set):
        """Genomic Prediction using XgBOOST with cross validation."""

        t1 = time.time()
        dataframes = []
        df_bim = pd.read_csv(bed_file.replace('bed', 'bim'), delimiter='\t', header=None)
        df_bim.columns = ['Chr', 'SNP', 'NA', 'ChrPos', 'NA', 'NA']

        for i in range(int(model_nr)):
            add_log('Model Iteration: ' + str(i + 1))

            # Split data into training and testing sets
            X_train, X_test, y_train, y_test = train_test_split(bed_fixed.read().val, pheno.read().val,
                                                                test_size=test_size)

            #scaler = StandardScaler()
            #X_train = scaler.fit_transform(X_train)
            xgb_model = xgb.XGBRegressor(n_estimators=estimators, learning_rate=0.1, max_depth=max_dep_set)
            xgb_model.fit(X_train, y_train)

            # for prediction, we need all genotypes
            bed_data = Bed(str(bed_file), count_A1=False, chrom_map=chrom_mapping)
            bed_data2 = bed_data.iid

            pheno_data2 = pheno.iid
            #predicted_values = xgb_model.predict(bed_fixed.read().val)
            predicted_values = xgb_model.predict(bed_data.read().val)

            df_gp = pd.DataFrame(bed_data2, columns=['ID1', 'BED_ID2'])
            df_gp['Predicted_Value'] = predicted_values

            df_pheno = pd.DataFrame(pheno_data2, columns=['ID1', 'Pheno_ID2'])
            df_pheno['Pheno_Value'] = pheno.read().val
            merged_df = df_gp.merge(df_pheno, on='ID1', how='outer')
            merged_df['Predicted_Value'] = merged_df['Predicted_Value'].astype(float)
            merged_df['Predicted_Value'] = merged_df['Predicted_Value'].round(5)
            dataframes.append(merged_df)


        df_all = self.helper.merge_gp_models(dataframes)
        df_all.to_csv(genomic_predict_name, index=False)

        t2 = time.time()
        t3 = round((t2 - t1) / 60, 2)
        add_log('Final run time (minutes): ' + str(t3))
        return df_all


    def start_gwas(self, bed_file, pheno_file, chrom_mapping, algorithm, add_log, test_size, model_nr, estimators, leave_chr_set,
                   max_dep_set, gwas_result_name, genomic_predict, genomic_predict_name):
        """Depreciated."""

        #print (chrom_mapping)
        #print (test_size, estimators, model_nr)
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

            elif algorithm == 'Linear regression':
                df = single_snp_linreg(test_snps=bed_fixed, pheno=pheno, output_file_name=gwas_result_name)


            elif algorithm == 'Random Forest (AI)':
                dataframes = []
                for i in range(int(model_nr)):
                    add_log('Model Iteration: ' + str(i +1))
                    df_bim = pd.read_csv(bed_file.replace('bed', 'bim'), delimiter='\t', header=None)
                    df_bim.columns = ['Chr', 'SNP', 'NA', 'ChrPos', 'NA', 'NA']
                    df = self.gwas_ai.run_random_forest(bed_fixed.read().val, pheno.read().val, df_bim, test_size,
                                                  estimators, gwas_result_name, bed_gp, pheno_gp, genomic_predict,
                                                  genomic_predict_name, model_nr)
                    #print (df)
                    df['Chr'] = df['Chr'].replace(chrom_mapping)
                    dataframes.append(df)
                # We merge the models and calculate the sum effect
                if genomic_predict:
                    df = df
                else:
                    df = self.helper.merge_models(dataframes)
                #df = df.sort_values(by='PValue', ascending=False)
            elif algorithm == 'XGBoost (AI)':
                dataframes = []
                for i in range(model_nr):
                    add_log('Model Iteration: ' + str(i + 1))
                    df_bim = pd.read_csv(bed_file.replace('bed', 'bim'), delimiter='\t', header=None)
                    df_bim.columns = ['Chr', 'SNP', 'NA', 'ChrPos', 'NA', 'NA']
                    #print(pheno_gp.iid[0:10], bed_gp.iid[0:10])
                    df = self.gwas_ai.run_xgboost(bed_fixed.read().val, pheno.read().val, df_bim, test_size,
                                                  estimators, str(i) + gwas_result_name, bed_gp, pheno_gp, genomic_predict,
                                                  genomic_predict_name, max_dep_set, model_nr)
                    #print (df)
                    #df['Chr'] = df['Chr'].replace(chrom_mapping)

                    dataframes.append(df)
                if genomic_predict:
                    df = df
                    #print (df)
                else:
                    df['Chr'] = df['Chr'].replace(chrom_mapping)
                    df = self.helper.merge_models(dataframes)
                #print (df)
                #df = df.sort_values(by='PValue', ascending=False)
                #print(df)

            elif algorithm == 'GP_LMM':
                #dataframes = []

                #for i in range(model_nr):
                    #print (i)
                df = pd.read_csv(bed_file.replace('bed', 'bim'), delimiter='\t', header=None)
                snp_ids = df.iloc[:, 1].tolist()
                #df = self.gwas_ai.run_xgboost(bed_fixed.read().val, pheno.read().val, snp_ids, test_size,
                                             # estimators, str(i) + gwas_result_name, bed_gp, pheno_gp, genomic_predict,
                                             # genomic_predict_name, max_dep_set, model_nr)
                df = self.gwas_ai.run_lmm_gp(bed_fixed, pheno, snp_ids, test_size,
                                              estimators,gwas_result_name, bed_gp, pheno_gp, genomic_predict,
                                              genomic_predict_name, max_dep_set, model_nr)
                    #dataframes.append(df)
                if genomic_predict:
                    df = df
                #else:
                    #df = self.helper.merge_models(dataframes)

            t2 = time.time()
            t3 = round((t2-t1)/ 60, 2)
            if algorithm == 'XGBoost (AI)' or algorithm == 'Random Forest (AI)' or algorithm == 'GP_LMM':
                if genomic_predict:
                    df = df.sort_values(by='Pheno_Value', ascending=False)
                    #df.dropna(subset=['Pheno_Value'], inplace=True)
                else:
                    df = df.sort_values(by='PValue', ascending=False)
                    #df.dropna(subset=['PValue'], inplace=True)


            #print (df)
            # we create one df for the plotting with ints as chr
            df_plot = df.copy(deep=True)

            if len(chrom_mapping) > 0:
                reversed_chrom_map = {value: key for key, value in chrom_mapping.items()}
                df["Chr"] = df["Chr"].apply(lambda x: reversed_chrom_map[x])

            df.to_csv(gwas_result_name, index=0)
            add_log('Final run time (minutes): ' + str(t3))
            return df, df_plot
        else:
            add_log(check_input_data[1], error=True)

    def plot_gwas(self, df, limit, algorithm, manhatten_plot_name, qq_plot_name, chrom_mapping):
        """Manhatten and qq-plot."""

        # Extract the top10 SNPs and use the value as significant marker label threshold
        if algorithm == 'FaST-LMM' or algorithm == 'Linear regression':

            # Remove rows where column 'A' has a value of 0.0

            df = df.sort_values(by=['Chr', 'ChrPos'])
            df = df[df['PValue'] != 0.0]

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

            flipped_dict = {value: key for key, value in chrom_mapping.items()}
            df['Chr'] = df['Chr'].astype(float).replace(flipped_dict)

            _ = gv.manhattanplot(data=df,chrom='Chr', pos="ChrPos", pv="PValue", snp="SNP", marker=".",color=['#4297d8', '#eec03c','#423496','#495227','#d50b6f','#e76519','#d580b7','#84d3ac'],
                              sign_marker_color="r", title="Manhatten Plot " + algorithm + '\n',
                              xlabel="Chromosome", ylabel=r"$-log_{10}{(P)}$", sign_line_cols=["#D62728", "#2CA02C"],
                              hline_kws={"linestyle": "--", "lw": 1.3},#, sign_marker_p=1e-9, is_annotate_topsnp=True,
                                 text_kws={"fontsize": 12, "arrowprops": dict(arrowstyle="-", color="k", alpha=0.6)},
                                 logp=True,ax=ax,xticklabel_kws={"rotation": "vertical"})

            plt.tight_layout(pad=1)
            plt.savefig(manhatten_plot_name)
            plt.savefig(manhatten_plot_name.replace('manhatten_plot', 'manhatten_plot_high'), dpi=1200)

            # Create QQ plot
            f, ax = plt.subplots(figsize=(6, 6), facecolor="w", edgecolor="k")
            #

            _ = gv.qqplot(data=df["PValue"], marker="o", title="QQ Plot " + algorithm + '\n',
                          xlabel=r"Expected $-log_{10}{(P)}$", ylabel=r"Observed $-log_{10}{(P)}$", ax=ax)
            # #ax.set(ylim=(0, 20), xlim=(0, 20))
            plt.tight_layout(pad=1)
            plt.savefig(qq_plot_name)
            plt.savefig(qq_plot_name.replace('qq_plot', 'qq_plot_high'), dpi=1200)


        else:
            plt_params = {
                "font.sans-serif": "Arial",
                "legend.fontsize": 10,
                "axes.titlesize": 14,
                "axes.labelsize": 12,
                "xtick.labelsize": 10,
                "ytick.labelsize": 10
            }
            plt.rcParams.update(plt_params)

            df = df.sort_values(by=['Chr', 'ChrPos'])
            df['Chr'] = df['Chr'].astype(int)
            df['ChrPos'] = df['ChrPos'].astype(int)

            flipped_dict = {value: key for key, value in chrom_mapping.items()}
            df['Chr'] = df['Chr'].astype(float).replace(flipped_dict)
            f, ax = plt.subplots(figsize=(12, 6), facecolor="w", edgecolor="k")
            algorithm2 = algorithm.replace(' (AI)', '')
            _ = gv.manhattanplot(data=df, chrom='Chr', pos="ChrPos", pv="PValue", snp="SNP", logp=False,
                                  title="Manhatten Plot " + algorithm2 + '\n',color=['#4297d8', '#eec03c','#423496','#495227','#d50b6f','#e76519','#d580b7','#84d3ac'],
                                  xlabel="Chromosome", ylabel=r"SNP effect", xticklabel_kws={"rotation": "vertical"})#, sign_marker_p=sign_limit, is_annotate_topsnp=True)
            plt.tight_layout(pad=1)
            plt.savefig(manhatten_plot_name)
            plt.savefig(manhatten_plot_name.replace('manhatten_plot', 'manhatten_plot_high'), dpi=1200)


    def plot_gp(self, df, gp_plot_name, algorithm):
        """Bland-Altman Plot for the real and predicted phenotype values."""


        # Calculate mean and difference (redundant here, but for demonstration)
        df['Mean'] = (df['Pheno_Value'] + df['Mean_Predicted_Value']) / 2
        df['Difference'] = df['Pheno_Value'] - df['Mean_Predicted_Value']

        # Plotting the Bland-Altman Plot
        plt.figure(figsize=(10, 6))
        plt.scatter(df['Mean'], df['Difference'], color='blue')
        #sns.scatterplot(x='Mean', y='Difference', data=df, color='blue')

        # Calculate and plot the mean difference
        mean_diff = df['Difference'].mean()
        plt.axhline(mean_diff, color='red', linestyle='--')

        # Calculate and plot the limits of agreement
        std_diff = df['Difference'].std()

        plt.axhline(mean_diff, color='red', linestyle='--', label='Mean Difference')
        plt.axhline(mean_diff + 1.96 * std_diff, color='green', linestyle='--', label='Upper Limit of Agreement')
        plt.axhline(mean_diff - 1.96 * std_diff, color='green', linestyle='--', label='Lower Limit of Agreement')
        plt.xlabel('Mean Value')
        plt.ylabel('Difference')
        algorithm2 = algorithm.replace(' (AI)', '')
        plt.title('Bland-Altman Plot (' + algorithm2 + ')')
        plt.legend()

        plt.tight_layout(pad=1)
        plt.savefig(gp_plot_name)
        plt.savefig(gp_plot_name.replace('Bland_Altman_plot', 'Bland_Altman_plot_high'), dpi=1200)


    def plot_gp_scatter(self, df, gp_plot_name_scatter, algorithm):
        """Regression Plot for the real and predicted phenotype values."""


        df = df.replace([np.inf, -np.inf], np.nan).dropna()

        # Pearson correlation
        corr, _ = pearsonr(df['Pheno_Value'], df['Mean_Predicted_Value'])
        corr_label = f"Pearson correlation: {corr:.2f}"

        # Plotting
        plt.figure(figsize=(10, 6))
        sns.scatterplot(data=df, x='Pheno_Value', y='Mean_Predicted_Value')
        sns.regplot(x='Pheno_Value', y='Mean_Predicted_Value', data=df, scatter=False, color='red')
        plt.title('Scatter Plot with Regression Line')
        plt.xlabel('Phenotype Value')
        plt.ylabel('Mean Predicted Value')
        #plt.text(5, max(df['Mean_Predicted_Value']) - 5, corr_label, fontsize=12, color='blue')
        algorithm2 = algorithm.replace(' (AI)', '')
        plt.title('Correlation Plot (' + algorithm2 + ')' + '\n' + corr_label)
        plt.tight_layout(pad=1)
        plt.savefig(gp_plot_name_scatter)
        plt.savefig(gp_plot_name_scatter.replace('GP_scatter_plot', 'GP_scatter_plot_high'), dpi=1200)
        #plt.show()


