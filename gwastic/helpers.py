import pandas as pd


class HELPERS:
    def duplicate_column(self, input_file, output_file):
        # Read the text file with one column
        df = pd.read_csv(input_file, header=None)
        df[1] = df[0]
        # Write the DataFrame with two columns to a new file
        df.to_csv(output_file, index=False, header=False, sep=' ')

    def replace_with_integers(self, input_file):
        print (input_file)
        mapping = {}  # Store the mapping of strings to integers
        current_integer = 1
        s = ''
        with open(input_file, 'r', errors="ignore") as infile:
            for line in infile:
                parts = line.strip().split('\t')
                col1_value = parts[0]

                # Check if the string in column 1 is already mapped to an integer
                if col1_value in mapping:
                    parts[0] = str(mapping[col1_value])
                else:
                    # If it's not mapped, assign the next integer and update the mapping
                    mapping[col1_value] = current_integer
                    parts[0] = str(current_integer)
                    current_integer += 1
                s += '\t'.join(parts) + '\n'
        with open(input_file, 'w', encoding='utf-8') as outfile:
            outfile.write(s)
        return mapping


