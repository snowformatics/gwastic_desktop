import pandas as pd
import numpy as np
from datetime import datetime
import os
import shutil
import configparser


class HELPERS:
    def duplicate_column(self, input_file, output_file):
        # Read the text file with one column
        df = pd.read_csv(input_file, header=None)
        df[1] = df[0]
        # Write the DataFrame with two columns to a new file
        df.to_csv(output_file, index=False, header=False, sep=' ')

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
                    if col1_value in mapping:
                        parts[0] = str(mapping[col1_value])
                    else:
                        # If it's not mapped, assign the next integer and update the mapping
                        mapping[col1_value] = current_integer
                        parts[0] = str(current_integer)
                        current_integer += 1
        return mapping

    def get_timestamp(self):
        """Get timestamp."""
        now = datetime.now()
        dt_string = now.strftime("%d%m%Y_%H%M%S")
        return dt_string

    def save_raw_data(self, bed, pheno):
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
                     qq_plot_name, algorithm, genomic_predict_name):
        """Stores the files of the analysis in the selected path."""
        ts = self.get_timestamp() + '_' + algorithm.replace(' ', '_')

        os.mkdir(os.path.join(save_dir, ts))
        save_dir = os.path.join(save_dir, ts)

        try:
            shutil.copyfile(os.path.join(current_dir, genomic_predict_name), os.path.join(save_dir, genomic_predict_name))

        except FileNotFoundError:
            shutil.copyfile(os.path.join(current_dir, manhatten_plot_name), os.path.join(save_dir, manhatten_plot_name))
            # We store also a trimmed version of single_snp with 10000 SNPs
            # todo if ai first sort
            df = pd.read_csv(gwas_result_name)
            first_10000_rows = df.head(10000)
            first_10000_rows.to_csv(gwas_result_name_top, index=False)

            shutil.copyfile(os.path.join(current_dir, gwas_result_name), os.path.join(save_dir, gwas_result_name))
            shutil.copyfile(os.path.join(current_dir, gwas_result_name_top), os.path.join(save_dir, gwas_result_name_top))

            if algorithm == "FaST-LMM:" or algorithm == "Linear regression":
                shutil.copyfile(os.path.join(current_dir, qq_plot_name), os.path.join(save_dir, qq_plot_name))
        return save_dir

    # def delete_files(self, current_dir, gwas_result_name, gwas_result_name_top, manhatten_plot_name, qq_plot_name,
    #                  algorithm):
    #     """Delete all temp files."""
    #     os.remove(os.path.join(current_dir, manhatten_plot_name))
    #     os.remove(os.path.join(current_dir, gwas_result_name))
    #     os.remove(os.path.join(current_dir, gwas_result_name_top))
    #     if algorithm == "FaST-LMM:" or algorithm == "Linear regression":
    #         os.remove(os.path.join(current_dir, qq_plot_name))








