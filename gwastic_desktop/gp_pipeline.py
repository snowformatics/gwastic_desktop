from gwastic_desktop.helpers import HELPERS
from gwastic_desktop.gwas_ai import GWASAI
from pysnptools.snpreader import Bed, Pheno
import time
from fastlmm.inference import FastLMM
from sklearn.model_selection import KFold
import pandas as pd
from sklearn.model_selection import train_test_split
import xgboost as xgb
from sklearn.ensemble import RandomForestRegressor
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import seaborn as sns

plt.switch_backend('Agg')


class GenomicPrediction:
    """GWAS class."""

    def __init__(self):
        self.gwas_ai = GWASAI()
        self.helper = HELPERS()

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
            # test = bed_fixed[test_indices, :]
            # Train the model
            fastlmm = FastLMM(GB_goal=2)
            fastlmm.fit(K0_train=train, y=pheno)
            # Test the model
            mean, covariance = fastlmm.predict(K0_whole_test=bed_data)
            df = pd.DataFrame({'Column1': mean.val[:, 0], 'Column2': mean.iid[:, 0]})

            # bed_data = bed_fixed.iid
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
        # df_bim = pd.read_csv(bed_file.replace('bed', 'bim'), delimiter='\t', header=None)
        # df_bim.columns = ['Chr', 'SNP', 'NA', 'ChrPos', 'NA', 'NA']

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

            # scaler = StandardScaler()
            # X_train = scaler.fit_transform(X_train)
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
        # df_bim = pd.read_csv(bed_file.replace('bed', 'bim'), delimiter='\t', header=None)
        # df_bim.columns = ['Chr', 'SNP', 'NA', 'ChrPos', 'NA', 'NA']
        # print (len(pheno.read().val))
        validation = False
        if validation:

            snp_data = bed_fixed.read().val
            pheno_data = pheno.read().val

            # Assuming x and y are your numpy arrays and they have the same length
            total_data_points = len(pheno_data)

            # Calculate the number of validation samples (5% of the total data)
            num_validation_samples = round(0.1 * total_data_points)

            # Generate random indices for the validation set
            indices = np.arange(total_data_points)
            np.random.shuffle(indices)
            val_indices = indices[:num_validation_samples]
            train_test_indices = indices[num_validation_samples:]

            # Create validation set
            snps_for_validation = snp_data[val_indices]
            pheno_for_validation = pheno_data[val_indices]
            ids_for_validation1 = bed_fixed.iid[val_indices]
            ids_for_validation2 = pheno.iid[val_indices]
            # print (ids_for_validation1)
            # print(ids_for_validation2, pheno_for_validation)

            # Remaining data for training/testing
            snps_for_training = snp_data[train_test_indices]
            pheno_for_training = pheno_data[train_test_indices]
            # print (bed_fixed.iid[train_test_indices])
            # np.savetxt('rest.csv', pheno_for_training, delimiter=',')
            # np.savetxt('val.csv', pheno_for_validation, delimiter=',')
        else:

            snps_for_training = bed_fixed.read().val
            pheno_for_training = pheno.read().val

        # print (len(snps_for_training), len(pheno_for_training))

        for i in range(int(model_nr)):
            add_log('Model Iteration: ' + str(i + 1))

            # Split data into training and testing sets
            # X_train, X_test, y_train, y_test = train_test_split(bed_fixed.read().val, pheno.read().val,
            # test_size=test_size)
            X_train, X_test, y_train, y_test = train_test_split(snps_for_training, pheno_for_training,
                                                                test_size=test_size)

            # scaler = StandardScaler()
            # X_train = scaler.fit_transform(X_train)
            xgb_model = xgb.XGBRegressor(n_estimators=estimators, learning_rate=0.1,
                                         max_depth=max_dep_set)  # , min_child_weight=3)
            xgb_model.fit(X_train, y_train)

            # for prediction, we need all genotypes
            bed_data = Bed(str(bed_file), count_A1=False, chrom_map=chrom_mapping)
            bed_data2 = bed_data.iid

            pheno_data2 = pheno.iid
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

        if validation:
            unique_ids = np.unique(ids_for_validation1.flatten())
            # Filter the DataFrame based on these IDs
            valdiation_df = df_all[df_all['ID1'].isin(unique_ids)]

            valdiation_df.to_csv(genomic_predict_name.replace('.csv', '_valdation.csv'), index=False)
            # print(filtered_df)
        else:
            valdiation_df = None

        t2 = time.time()
        t3 = round((t2 - t1) / 60, 2)
        add_log('Final run time (minutes): ' + str(t3))
        return df_all, valdiation_df

    def plot_gp(self, df, gp_plot_name, algorithm):
        """Bland-Altman Plot for the real and predicted phenotype values."""

        # Calculate mean and difference (redundant here, but for demonstration)
        df['Mean'] = (df['Pheno_Value'] + df['Mean_Predicted_Value']) / 2
        df['Difference'] = df['Pheno_Value'] - df['Mean_Predicted_Value']

        # Plotting the Bland-Altman Plot
        plt.figure(figsize=(10, 6))
        plt.scatter(df['Mean'], df['Difference'], color='blue')
        # sns.scatterplot(x='Mean', y='Difference', data=df, color='blue')

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
        plt.savefig(gp_plot_name.replace('Bland_Altman_plot', 'Bland_Altman_plot_high'), dpi=300)

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
        # plt.text(5, max(df['Mean_Predicted_Value']) - 5, corr_label, fontsize=12, color='blue')
        algorithm2 = algorithm.replace(' (AI)', '')
        plt.title('Correlation Plot (' + algorithm2 + ')' + '\n' + corr_label)
        plt.tight_layout(pad=1)
        plt.savefig(gp_plot_name_scatter)
        plt.savefig(gp_plot_name_scatter.replace('GP_scatter_plot', 'GP_scatter_plot_high'), dpi=300)
        # plt.show()
