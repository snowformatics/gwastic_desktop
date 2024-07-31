import pandas as pd

def create_val_sets():
    l = ["WCC_for_GWAS_Bgt96224_SL_blues5", "WCC_for_GWAS_FAL92315_SL_blues5", "WCC_for_GWAS_mean_PM1_PM2_SL_blues5"]
    path = "Z:/wheat/genbank20_GP/"
    nr = 5
    for i in l:
        # Step 1: Load the CSV data into a DataFrame
        df = pd.read_csv(path + i + ".txt", delimiter=' ', names=['ID', 'value'])

        for n in range(nr):
            # Step 2: Randomly select 500 rows
            random_rows = df.sample(n=100)

            # Step 3: Remove the selected rows from the original DataFrame
            remaining_df = df.drop(random_rows.index)

            remaining_df['ID2'] = remaining_df['ID']
            remaining_df = remaining_df[['ID', 'ID2', 'value']]
            #print (remaining_df)
            # Step 4: Save the two DataFrames to two separate CSV files
            random_rows.to_csv(path + i + '_val_' + str(n) + '.txt', sep=' ',index=False, header=None)
            remaining_df.to_csv(path + i + '_remaining_' + str(n) + '.txt', sep=' ',index=False, header=None)

def merge_val_sets():
    pass