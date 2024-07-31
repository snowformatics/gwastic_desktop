import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import uniform, randint
import matplotlib.pyplot as plt
import pandas as pd
import os
import glob

import geneview as gv


def plot_fi_all():

    # Path to the directory containing your CSV files
    #directory_path = "U:/Paper/gwastic/validation/validation_sets/barley/03012024_142206_Random_Forest_(AI)/iteration/"
    directory_path = "U:/Paper/gwastic/validation/validation_sets/barley/03012024_142719_XGBoost_(AI)/iteration/"
    # Replace with your actual directory path

    # Creating a pattern to match all CSV files in the directory
    file_pattern = os.path.join(directory_path, "*.csv")
    #print (file_pattern)
    # Finding all files in the directory that match the pattern
    csv_files = glob.glob(file_pattern)
    #print (csv_files)

    # Reading each CSV file into a DataFrame and storing them in a list
    dataframes = [pd.read_csv(file, names=['SNP', 'PValue'], sep=' ') for file in csv_files]
    #dataframes = [pd.read_csv(file,  sep=',') for file in csv_files]

    # Merging all dataframes
    merged_df = pd.concat(dataframes)
    #merged_df.columns = ['SNP', 'PValue']
    # Display the first few rows of the merged dataframe
    print (merged_df)
    filtered_df = merged_df[merged_df['PValue'] > 0]
    print (filtered_df)
    #filtered_df = merged_df
    filtered_df[['Chr', 'ChrPos']] = filtered_df['SNP'].str.split(':', expand=True)
    filtered_df = filtered_df.sort_values(by=['Chr', 'ChrPos'])
    #print(filtered_df.head())

    filtered_df['Chr'] = filtered_df['Chr'].astype(int)
    # chr_names = df['Chr'].unique()
    filtered_df['ChrPos'] = filtered_df['ChrPos'].astype(int)
    filtered_df.to_csv('all.csv', index=False)

    f, ax = plt.subplots(figsize=(12, 5), facecolor="w", edgecolor="k")

    _ = gv.manhattanplot(data=filtered_df, chrom='Chr', pos="ChrPos", pv="PValue", snp="SNP", logp=False,
                         title="GWAS Manhatten Plot XGB 10x" + '\n', color=['#3B5488', '#d5536f'],
                         xlabel="Chromosome", ylabel=r"Feature Importance")

    plt.tight_layout(pad=1)
    plt.savefig("mh_test.png", dpi=200)

    # # Filtering out the values that are above 0
    # filtered_df = merged_df[merged_df['value'] > 0]
    #
    # # Plotting the values
    # plt.figure(figsize=(10, 6))
    # plt.scatter(filtered_df['id'], filtered_df['value'])
    # plt.xlabel('ID')
    # plt.ylabel('Value')
    # plt.title('Values Above 0 in Merged Datasets')
    # plt.xticks(rotation=90)
    # plt.show()

#plot_fi_all()
# Create the data
# data = {
#     'ID': [6961, 6917, 6898, 6963, 6913, 6904, 6908, 6970, 6971, 7521, 7525],
#     'Predicted_Value': [0.95, 0.99, 1.00, 1.00, 1.00, 0.99995726, 0.99994034, 0.99981785, 0.9473342, 0.68087214, 0.83025223],
#     'Pheno_Value': [1, 1, 1, 1, 1, None, None, None, None, None, None],
#     'difference': [0.05, 0.01, 0.00, 0.00, 0.00, None, None, None, None, None, None]
# }
#
# df = pd.DataFrame(data)

# df = pd.read_csv('test.txt', delimiter='\t')
#
#
# # Create the scatter plot
# fig, ax = plt.subplots()
#
# # Plot the data points for Pheno_Value
# ax.scatter(df['ID'], df['Pheno_Value'], color='blue', label='Pheno_Value')
#
# # Only plot the Predicted_Value where the Pheno_Value is NaN
# mask = df['Pheno_Value'].isna()
# ax.scatter(df['ID'][mask], df['Predicted_Value'][mask], color='yellow', label='Predicted_Value')
#
# # Add the error bars only when difference is not 0
# non_zero_difference_mask = df['difference'] != 0
# ax.errorbar(df['ID'][non_zero_difference_mask], df['Predicted_Value'][non_zero_difference_mask], yerr=df['difference'][non_zero_difference_mask], fmt='none', color='red', capsize=5)
#
# # Add labels and legend
# ax.set_xlabel('ID')
# ax.set_ylabel('Value')
# ax.legend()
#
# plt.show()


# print (df)
# # Create the scatter plot
# fig, ax = plt.subplots()
#
# # Plot the data points
# ax.scatter(df['ID'], df['Predicted_Value'], color='yellow', label='Predicted_Value')
# ax.scatter(df['ID'], df['Pheno_Value'], color='blue', label='Pheno_Value')
#
# # Add the error bars
# ax.errorbar(df['ID'], df['Predicted_Value'], yerr=df['difference'], fmt='none', color='red', capsize=5)
#
# # Add labels and legend
# ax.set_xlabel('ID')
# ax.set_ylabel('Value')
# ax.legend()
#
# plt.show()


# import pandas as pd
# import matplotlib.pyplot as plt
# import seaborn as sns
#
# # Sample dataframe
# data = {
#     'Phenotype': ['A', 'B', 'A', 'C', 'B', 'A', 'C', 'B', 'B', 'C', 'A'],
# }
# df = pd.DataFrame(data)
#
# # Create a countplot (bar plot of value counts)
# plt.figure(figsize=(10, 6))
# sns.countplot(x='Phenotype', data=df, order=df['Phenotype'].value_counts().index)
# plt.title('Phenotype Distribution')
# plt.ylabel('Count')
# plt.xlabel('Phenotype')
# plt.show()


def sns_plot():
    # Simulate DataFrame
    # df = pd.DataFrame({
    # 'rsid'  : ['rs{}'.format(i) for i in np.arange(10000)],
    # 'chrom' : [i for i in randint.rvs(1,23+1,size=10000)],
    # 'pos'   : [i for i in randint.rvs(0,10**5,size=10000)],
    # 'pval'  : uniform.rvs(size=10000)})
    # df['-logp'] = -np.log10(df.pval)
    # df = df.sort_values(['chrom','pos'])

    df = pd.read_csv("single_snp.csv", delimiter='\t')
    df['-logp'] = -np.log10(df['PValue'])
    df = df.sort_values(['Chr','ChrPos'])

    df.reset_index(inplace=True, drop=True)
    df['i'] = df.index

    # Generate Manhattan plot: (#optional tweaks for relplot: linewidth=0, s=9)
    plot = sns.relplot(data=df, x='i', y='-logp', aspect=3.7,
                       hue='Chr', palette = 'bright', legend=None)
    chrom_df=df.groupby('Chr')['i'].median()
    plot.ax.set_xlabel('Chr')
    plot.ax.set_xticks(chrom_df)
    plot.ax.set_xticklabels(chrom_df.index)
    plot.fig.suptitle('Manhattan plot')
    plt.show()

# import matplotlib.pyplot as plt
# import geneview as gv
#
# df = pd.read_csv("gwas_results.csv", delimiter=' ')
# df.columns = ['SNP', 'PValue']
# df[['Chr', 'ChrPos']] = df['SNP'].str.split(':', expand=True)
#
#
# #df = df.sort_values(by=['Chr', 'ChrPos'])
# df['Chr'] = df['Chr'].astype(int)
# # chr_names = df['Chr'].unique()
# df['ChrPos'] = df['ChrPos'].astype(int)
# ax = gv.manhattanplot(data=df, chrom='Chr', pos="ChrPos", pv="PValue", snp="SNP", logp=False, title="GWAS Manhatten Plot",
#                       xlabel="Chromosome", ylabel=r"Feature Importance")
# plt.show()
# # df = df.head(limit)
#
# # common parameters for plotting
# plt_params = {
#     "font.sans-serif": "Arial",
#     "legend.fontsize": 14,
#     "axes.titlesize": 18,
#     "axes.labelsize": 16,
#     "xtick.labelsize": 14,
#     "ytick.labelsize": 14
# }
# plt.rcParams.update(plt_params)
#
# # Create a manhattan plot
# f, ax = plt.subplots(figsize=(12, 5), facecolor="w", edgecolor="k")
# # xtick = set(["chr" + i for i in list(map(str, chr_names))])
# _ = gv.manhattanplot(data=df, chrom='Chr', pos="ChrPos", pv="PValue", snp="SNP", marker="."
#                      , logp=False,
#                      sign_marker_color="r", title="GWAS Manhatten Plot " + '\n',
#                      # xtick_label_set=xtick,
#                      xlabel="Chromosome", ylabel=r"Feature Importance", sign_line_cols=["#D62728", "#2CA02C"],
#                      #hline_kws={"linestyle": "--", "lw": 1.3},
#                      text_kws={"fontsize": 12, "arrowprops": dict(arrowstyle="-", color="k", alpha=0.6)},
#                      ax=ax)
# plt.tight_layout(pad=1)
# plt.savefig('test.png', dpi=100)
#
# df = df.sort_values(by=['PValue'])
# print (df)
#
#



# set parameter show=True, if you want view the image instead of saving

#sns_plot()
#
# # import libraries
# from pandas import DataFrame
# from scipy.stats import uniform
# from scipy.stats import randint
# import numpy as np
# import matplotlib.pyplot as plt
# #import seaborn as sns
# #sns.set_theme()
# # sample data
# # df = DataFrame({'gene' : ['gene-%i' % i for i in np.arange(1000000)],
# # 'pvalue' : uniform.rvs(size=1000000),
# # 'Chr' : ['ch-%i' % i for i in randint.rvs(0,12,size=1000000)]})
# # # -log_10(pvalue)
# # df['-logp'] = -np.log10(df.pvalue)
# # print (df)
# df = pd.read_csv("single_snp.csv", delimiter='\t')
# df['-logp'] = -np.log10(df['PValue'])
# df['Chr'] = df['Chr'].apply(lambda x: f"chr{int(x)}")
# print (df)
# #df = df.sort_values(['Chr', 'ChrPos'])
#
#
# df['Chr'] = df['Chr'].astype('category')
#
# #df['Chr']= df['Chr'].cat.set_categories(['ch-%i' % i for i in range(12)], ordered=True)
# df['Chr']= df['Chr'].cat.set_categories(['chr%i' % i for i in range(12)], ordered=True)
#
# #print (df)
# df = df.sort_values('Chr')
#
#
# # How to plot gene vs. -log10(pvalue) and colour it by chromosome?
# df['ind'] = range(len(df))
# df_grouped = df.groupby(('Chr'))
# print (df_grouped)
# # manhattan plot
# fig = plt.figure(figsize=(14, 8)) # Set the figure size
# ax = fig.add_subplot(111)
# colors = ['#FFD1DC','#A2CFFE','#B0E57C', '#FFF79A']
# x_labels = []
# x_labels_pos = []
# #print (df)
# #df = df.reset_index()
# for num, (name, group) in enumerate(df_grouped):
#     group.plot(kind='scatter', x='ind', y='-logp',color=colors[num % len(colors)], ax=ax)
#     x_labels.append(name)
#     #print (group['ind'])
#     #print (group['ind'].iloc[-1])
#     #x_labels_pos.append((group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0])/2))
# #ax.set_xticks(x_labels_pos)
# ax.set_xticklabels(x_labels)
#
# # set axis limits
# ax.set_xlim([0, len(df)])
# ax.set_ylim([0, 20])
#
# # x axis label
# ax.set_xlabel('Chromosome')
#
# # show the graph
# plt.show()

import matplotlib.pyplot as plt
import geneview as gv

# load data
df = gv.load_dataset("gwas")
print (df)
# Plot a basic manhattan plot with horizontal xtick labels and the figure will display in screen.
ax = gv.manhattanplot(data=df)
plt.show()


