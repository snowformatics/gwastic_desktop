import pandas as pd

class HELPERS:

    def duplicate_column(self, input_file, output_file):
        # Read the text file with one column
        df = pd.read_csv(input_file, header=None)

        # Duplicate the column
        df[1] = df[0]

        # Write the DataFrame with two columns to a new file
        df.to_csv(output_file, index=False, header=False, sep=' ')

# Replace 'input.txt' with the path to your input file
# Replace 'output.txt' with the desired name for the output file
