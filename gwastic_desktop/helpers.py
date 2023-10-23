import pandas as pd
from datetime import datetime
import os
import shutil

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
        now = datetime.now()
        dt_string = now.strftime("%d%m%Y_%H%M%S")
        print (dt_string)
        return dt_string

    def save_results(self, current_dir, save_dir):
        ts = self.get_timestamp()

        # try:
        #     os.mkdir(os.path.join(save_dir, ts))
        #     save_dir = os.path.join(save_dir, ts)
        # except:
        #     save_dir = save_dir
        # #
        os.mkdir(os.path.join(save_dir, ts))
        save_dir = os.path.join(save_dir, ts)
        shutil.copyfile(os.path.join(current_dir, "manhatten.png"), os.path.join(save_dir, "manhatten.png"))
        shutil.copyfile(os.path.join(current_dir, "qq.png"), os.path.join(save_dir, "qq.png"))
        # We store also a trimmed version of single_snp with 10000 SNPs
        df = pd.read_csv("single_snp.csv")
        first_10000_rows = df.head(10000)
        first_10000_rows.to_csv("single_snp_top10000.csv", index=False)
        shutil.copyfile(os.path.join(current_dir, "single_snp.csv"), os.path.join(save_dir, "single_snp.csv"))
        shutil.copyfile(os.path.join(current_dir, "single_snp_top10000.csv"), os.path.join(save_dir, "single_snp_top10000.csv"))







