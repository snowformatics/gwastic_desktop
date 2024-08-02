import pandas as pd
import numpy as np
from datetime import datetime
import os
import shutil
import configparser


class HELPERS:
    # def duplicate_column(self, input_file, output_file):
    #     """Depreciated."""
    #     # Read the text file with one column
    #     df = pd.read_csv(input_file, header=None)
    #     df[1] = df[0]
    #     # Write the DataFrame with two columns to a new file
    #     df.to_csv(output_file, index=False, header=False, sep=' ')

    def replace_with_integers(self, input_file):
        """Replace string chromosome names with integers."""
        mapping = {}  # Store the mapping of strings to integers
        current_integer = 1.0
        with open(input_file, 'r', errors="ignore") as infile:
            for line in infile:
                parts = line.strip().split('\t')
                col1_value = parts[0]

                try:
                    col1_value = int(col1_value)
                except ValueError:
                    # Check if the string in column 1 is already mapped to an integer
                    #print(col1_value)
                    if col1_value in mapping:
                        parts[0] = str(mapping[col1_value])
                    else:
                        # If it's not mapped, assign the next integer and update the mapping
                        mapping[col1_value] = current_integer
                        parts[0] = str(current_integer)
                        current_integer += 1
        #print (mapping)
        return mapping

    def get_timestamp(self):
        """Get timestamp."""
        now = datetime.now()
        dt_string = now.strftime("%d%m%Y_%H%M%S")
        return dt_string

    def save_raw_data(self, bed, pheno):
        """Save the SNP matrix."""
        np.save('snp', bed.read().val)
        np.savez_compressed('snp.npz', bed.read().val)
        np.save('pheno', pheno.read().val)

    def save_settings(self, default_path):
        """Save the user settings."""
        config = configparser.ConfigParser()
        # Add sections and settings
        config['DefaultSettings'] = {'path': default_path, 'algorithm':'FaST-LMM'}
        # Write to a file
        with open('settings.ini', 'w') as configfile:
            config.write(configfile)

    def get_settings(self, setting):
        """Get the user settings."""
        config = configparser.ConfigParser()
        config.read('settings.ini')
        # Accessing values from DefaultSettings
        default_path= config['DefaultSettings'][setting]
        return default_path

    def save_results(self, current_dir, save_dir, gwas_result_name, gwas_result_name_top, manhatten_plot_name,
                     qq_plot_name, algorithm, genomic_predict_name, gp_plot_name, gp_plot_name_scatter, add_log,
                     settings_lst, pheno_stats_name, geno_stats_name ):
        """Stores the files of the analysis in the selected path."""
        ts = self.get_timestamp() + '_' + algorithm.replace(' ', '_')
        save_dir = os.path.join(save_dir, ts)

        try:
            os.mkdir(save_dir)
        except OSError:
            add_log('Can not create folder. Please select a valid directory.', error=True)
        try:
            # Create a gwas output file with top 10000 SNPs
            df = pd.read_csv(gwas_result_name)
            first_10000_rows = df.head(10000)
            first_10000_rows.to_csv(gwas_result_name_top, index=False)
        except FileNotFoundError:
            pass
        # Copy all results files
        src_files = [gwas_result_name, gwas_result_name_top, genomic_predict_name,
                     #manhatten_plot_name, qq_plot_name, gp_plot_name, gp_plot_name_scatter,
                     manhatten_plot_name.replace('manhatten_plot', 'manhatten_plot_high'),
                     qq_plot_name.replace('qq_plot', 'qq_plot_high'),
                     gp_plot_name_scatter.replace('GP_scatter_plot', 'GP_scatter_plot_high'),
                     gp_plot_name.replace('Bland_Altman_plot', 'Bland_Altman_plot_high'),
                     genomic_predict_name.replace('.csv', '_valdation.csv'),
                     pheno_stats_name, geno_stats_name]
        for src_file in src_files:
            if os.path.exists(os.path.join(current_dir, src_file)):
                shutil.copy(os.path.join(current_dir,src_file), os.path.join(save_dir,src_file))
                add_log(f"File saved: {src_file}")
            else:
                pass
        # Create a log file with all settings
        log_file = open(save_dir + '/log.txt', 'w')
        log_file.write('Algorithm: ' + settings_lst[0])
        log_file.write('\nBed file used: ' + settings_lst[1])
        log_file.write('\nPheno file used: ' + settings_lst[2])
        log_file.write('\nTraining size: ' + str(100-(100*float(settings_lst[3]))))
        log_file.write('\nNr of trees: ' + str(settings_lst[4]))
        log_file.write('\nNr of models: ' + str(settings_lst[5]))
        log_file.write('\nMax depth: ' + str(settings_lst[6]))
        return save_dir

    def merge_gp_models(self, dataframes):
        """Merge all GP ML models."""

        # """Combine RF or XG models and calculate the mean and SD of the Predicted_Value."""
        # df_combined = pd.concat(dataframes)
        #
        # # Grouping by 'ID1' and 'BED_ID2' to calculate mean and SD of 'Predicted_Value'
        # df_grouped = df_combined.groupby(['ID1', 'BED_ID2'])['Predicted_Value'].agg(['mean', 'std']).reset_index()
        # df_grouped = df_grouped.rename(columns={'mean': 'Mean_Predicted_Value', 'std': 'SD_Predicted_Value'})
        #
        # # Merging the grouped results with the first dataframe in the list
        # df_result = pd.merge(dataframes[0], df_grouped, on='ID1', how='left')
        # df_result.to_csv(('out.txt'))
        # df_result = df_result.drop(['Predicted_Value_x', 'BED_ID2_y', 'Pheno_ID2'], axis=1)
        # df_result = df_result.rename(columns={'Predicted_Value_y': 'Mean_Predicted_Value'})
        #
        # # Calculate the absolute difference between 'Pheno_Value' and 'Mean_Predicted_Value'
        # df_result['Difference'] = (df_result['Pheno_Value'] - df_result['Mean_Predicted_Value']).abs()
        #
        # # Rounding the values
        # df_result['Difference'] = df_result['Difference'].round(decimals=3)
        # df_result['Mean_Predicted_Value'] = df_result['Mean_Predicted_Value'].round(decimals=3)
        # df_result['SD_Predicted_Value'] = df_result['SD_Predicted_Value'].round(decimals=3)
        #
        # return df_result

        df_combined = pd.concat(dataframes)
        df_result = df_combined.groupby(['ID1', 'BED_ID2'])['Predicted_Value'].mean().reset_index()
        df_result = pd.merge(dataframes[0], df_result, on='ID1', how='left')
        df_result = df_result.drop(['Predicted_Value_x', 'BED_ID2_y', 'Pheno_ID2'], axis=1)
        df_result = df_result.rename(columns={'Predicted_Value_y': 'Mean_Predicted_Value'})
        df_result['Difference'] = (df_result['Pheno_Value'] - df_result['Mean_Predicted_Value']).abs()
        df_result['Difference'] = df_result['Difference'].round(decimals=3)
        df_result['Mean_Predicted_Value'] = df_result['Mean_Predicted_Value'].round(decimals=3)
        return df_result

    def merge_models(self, dataframes, method):
        """Merge all GWAS ML models."""

        # """Combine RF or XG models and calculate the sum of the SNP effect."""
        # df_combined = pd.concat(dataframes)
        # df_combined = df_combined[df_combined['PValue'] > 0]
        # # Grouping by 'snp' and summing the values
        # df2 = df_combined.groupby('SNP')['PValue'].sum().reset_index()
        # df_result_sum = pd.merge(df_combined, df2, on='SNP', how='left')
        # df_result_sum = df_result_sum.drop(['PValue_x'], axis=1).drop_duplicates()
        # df_result_sum = df_result_sum.rename(columns={'PValue_y': 'PValue'})
        # df_result_sum['Chr'] = df_result_sum['Chr'].astype(int)
        # df_result_sum['ChrPos'] = df_result_sum['ChrPos'].astype(int)
        # df_result_sum = df_result_sum.sort_values(by=['Chr', 'ChrPos'])
        # return df_result_sum

        """Combine RF or XG models and calculate the sum and SD of the SNP effect."""
        df_combined = pd.concat(dataframes)

        df_combined = df_combined[df_combined['PValue'] > 0]
        #df_combined.to_csv('com.txt')
        # Grouping by 'SNP' and calculating the sum and SD of PValues
        df_grouped = df_combined.groupby(['SNP', 'Chr', 'ChrPos'])['PValue'].agg([method, 'std']).reset_index()
        #df_grouped.to_csv('group.txt')
        df_grouped = df_grouped.rename(columns={method: 'PValue_x', 'std': 'PValue_sd'})
        #df_result = pd.merge(df_grouped, df_combined, on='SNP', how='left')
        #df_result.to_csv('res1.txt')
        #df_result = df_result.drop(['PValue_x'], axis=1).drop_duplicates()
        #df_result.to_csv('res2.txt')

        df_result = df_grouped
        df_result = df_result.rename(columns={'PValue_x': 'PValue'})
        df_result['Chr'] = df_result['Chr'].astype(int)
        df_result['ChrPos'] = df_result['ChrPos'].astype(int)
        df_result = df_result.sort_values(by=['Chr', 'ChrPos'])
        #df_result.to_csv('res1.txt')

        return df_result











