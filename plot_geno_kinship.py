import numpy as np
from fastlmm.util import example_file # Download and return local file name
from pysnptools.snpreader import Bed, Pheno
import pysnptools.util as pstutil
from pysnptools.kernelreader import SnpKernel
from pysnptools.standardizer import Unit
import seaborn as sns
import matplotlib.pyplot as plt


bed_fn = "//filer-5/agruppen/BIT/lueck/gwas_test_data/gwastic_paper/small_set/example.bed"
pheno_fn = "pheno_gwas.csv"
bed_fn = "//filer-5/agruppen/BIT/lueck/gwas_test_data/gwastic_paper/WGS300_005_0020.bed"
pheno_fn = "pheno_gwas2.csv"


def filter_out_missing(bed):
    """Filter out missing."""

    sid_batch_size = 1000  # number of SNPs to read at a time, decrease in case we run oput of memory
    # list of all NaN columns
    all_nan = []

    for sid_start in range(0, bed.sid_count, sid_batch_size):

        # read a batch of SNPs
        snp_data_batch = bed[:, sid_start:sid_start + sid_batch_size].read()
        # find the columns that are all NaN (missing)
        nan_in_batch = np.isnan(snp_data_batch.val).all(axis=0)
        all_nan.append(nan_in_batch)

    # concatenate the list of all NaN columns into a single boolean array
    all_nan = np.concatenate(all_nan, axis=0)
    # print the number of all NaN columns

    # remove these SNPs from consideration
    bed_fixed = bed[:, ~all_nan]
    return bed_fixed


bed = Bed(bed_fn, count_A1=True).read()
pheno = Pheno(pheno_fn).read()
print (bed.shape, pheno.shape, type(bed), type(pheno))
bed, pheno = pstutil.intersect_apply([bed, pheno])
bed = filter_out_missing(bed)
bed = bed.read()
pheno = pheno.read()

print (bed.shape, pheno.shape, type(bed), type(pheno))


#print (pheno.val)
#bed = filter_out_missing(bed)
#print (bed2.val)
print (bed.shape, pheno.shape, type(bed))


print(f"Bed dimensions: {bed.iid_count} x {bed.sid_count}")

# names of the samples

print(f"Bed dimensions: {bed.iid_count} x {bed.sid_count}")
# sample bed values
#print(f"sample bed values: \n{bed.val[0:5,0:5]}")
#print(f"sample bed values: \n{bed.val()}")

# standardize the bed data
# that is: Make it mean 0 and sd 1 with any missing values filled in with 0
bed_standardized = bed.standardize()
#print(f"sample bed values: \n{bed.val[0:5,0:5]}")

# create kinship matrix
kernel = SnpKernel(bed_standardized, standardizer=Unit())
kernel = kernel.read().standardize()  # defaults to DiagKtoN standardize
print(f"sample kernel values: \n{kernel.val[0:5,0:5]}")

# Behind the scenes this is the same as:
# kernel2 = snps * snps.T / snps.shape[1]
kernel2 = bed_standardized.val @ bed_standardized.val.T / bed_standardized.sid_count
print(f"sample kernel2 values: \n{kernel2[0:5,0:5]}")


# # Extract phenotype values
# y = pheno.val.flatten()
#
# # Add a column of ones to the covariates to include an intercept
# X = np.ones((len(y), 1))
#
# # Fit the mixed linear model
# model = sm.MixedLM(endog=y, exog=X, groups=np.arange(len(y)), exog_re=kernel2)
# result = model.fit()
#
# # Extract variance components
# vc_random = result.cov_re.iloc[0, 0]
# vc_residual = result.scale
#
# # Calculate heritability
# heritability = vc_random / (vc_random + vc_residual)
# heritability_percent = heritability * 100
#
# print(f"Heritability: {heritability:.4f}")
# print(f"Heritability (%): {heritability_percent:.2f}%")

# plot the kinship matrix
import pylab
# Plotting the kinship matrix
# Extract sample names
sample_names = [''.join(iid) for iid in bed.iid[:, 1]]

# Create a DataFrame for seaborn
import pandas as pd
kinship_df = pd.DataFrame(kernel2, index=sample_names, columns=sample_names)

# Plotting the hierarchically-clustered heatmap with adjusted size and labels
g = sns.clustermap(kinship_df, cmap='viridis', figsize=(15, 15), xticklabels=1, yticklabels=1)

# Adjust the title positioning
plt.subplots_adjust(top=0.95)
g.fig.suptitle('Hierarchically-Clustered Kinship Matrix', fontsize=16)

# Adjust the font size for labels
plt.setp(g.ax_heatmap.get_xticklabels(), fontsize=10, rotation=90)
plt.setp(g.ax_heatmap.get_yticklabels(), fontsize=10, rotation=0)

# Move the color bar to the side
#g.cax.set_position([0.95, 0.2, 0.03, 0.45])
plt.show()
