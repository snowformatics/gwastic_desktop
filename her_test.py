import numpy as np
from pysnptools.snpreader import Bed, Pheno
import pysnptools.util as pstutil
from pysnptools.standardizer import Unit
from pysnptools.kernelreader import SnpKernel
import matplotlib.pyplot as plt
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

def h1():
    # File paths
    bed_fn = "//filer-5/agruppen/BIT/lueck/gwas_test_data/gwastic_paper/small_set/example.bed"
    pheno_fn = "pheno_gwas.csv"
    bed_fn = "//filer-5/agruppen/BIT/lueck/gwas_test_data/gwastic_paper/WGS300_005_0020.bed"
    #pheno_fn = "pheno_gwas2.csv"


    # Read BED file
    bed = Bed(bed_fn, count_A1=True).read()
    print(f"Bed dimensions: {bed.iid_count} x {bed.sid_count}")

    # Read phenotype file
    pheno = Pheno(pheno_fn).read()

    # Intersect bed and phenotype data
    bed, pheno = pstutil.intersect_apply([bed, pheno])
    bed = filter_out_missing(bed)
    bed = bed.read()
    pheno = pheno.read()
    print(f"Bed dimensions: {bed.iid_count} x {bed.sid_count}")

    # Standardize the bed data
    bed_standardized = bed.standardize(Unit())
    print(f"Sample bed values: \n{bed_standardized.val[0:5,0:5]}")

    # Create kinship matrix
    kernel = SnpKernel(bed_standardized, standardizer=Unit())
    kernel = kernel.read().standardize()  # defaults to DiagKtoN standardize
    K = kernel.val
    print(f"Sample kernel values: \n{K[0:5,0:5]}")

    # Extract phenotype values
    y = pheno.val.flatten()

    # Calculate genetic variance (Vg) and residual variance (Ve)
    Vg = np.var(K @ y)
    Ve = np.var(y - K @ y)

    # Calculate heritability (h2)
    heritability = Vg / (Vg + Ve)
    heritability_percent = heritability * 100

    print(f"Heritability: {heritability:.4f}")
    print(f"Heritability (%): {heritability_percent:.2f}%")

def h2():
    import numpy as np
    from pysnptools.snpreader import Bed, Pheno
    import pysnptools.util as pstutil
    from pysnptools.standardizer import Unit
    from pysnptools.kernelreader import SnpKernel
    from sklearn.decomposition import PCA

    # File paths
    bed_fn = "//filer-5/agruppen/BIT/lueck/gwas_test_data/gwastic_paper/small_set/example.bed"
    pheno_fn = "pheno_gwas.csv"
    bed_fn = "//filer-5/agruppen/BIT/lueck/gwas_test_data/gwastic_paper/WGS300_005_0020.bed"
    #pheno_fn = "pheno_gwas2.csv"
    pheno_fn = "//filer-5/agruppen/BIT/lueck/NH_Paper/barley/phenoytpes/random_blues.txt"

    # Read BED file
    bed = Bed(bed_fn, count_A1=True).read()
    #print(f"Bed dimensions: {bed.iid_count} x {bed.sid_count}")

    # Read phenotype file
    pheno = Pheno(pheno_fn).read()

    # Intersect bed and phenotype data
    bed, pheno = pstutil.intersect_apply([bed, pheno])
    bed = filter_out_missing(bed)
    bed = bed.read()
    pheno = pheno.read()

    print(f"Bed dimensions: {bed.iid_count} x {bed.sid_count}")
    #print (pheno.val)
    # Standardize the bed data
    bed_standardized = bed.standardize(Unit())
    #print(f"Sample bed values: \n{bed_standardized.val[0:5,0:5]}")

    # Create kinship matrix
    kernel = SnpKernel(bed_standardized, standardizer=Unit())
    kernel = kernel.read().standardize()  # defaults to DiagKtoN standardize
    K = kernel.val
    #print(f"Sample kernel values: \n{K[0:5,0:5]}")

    # Extract genotype and phenotype values
    X = bed_standardized.val
    y = pheno.val.flatten()

    # Calculate principal components
    pca = PCA(n_components=10)
    pcs = pca.fit_transform(X)
    print(f"Principal Components shape: {pcs.shape}")

    # Add principal components as covariates
    X_with_pcs = np.hstack((pcs, np.ones((pcs.shape[0], 1))))  # Include an intercept
    print(f"X_with_pcs shape: {X_with_pcs.shape}")

    # Adjust phenotype for population structure
    y_adjusted = y - X_with_pcs @ np.linalg.lstsq(X_with_pcs, y, rcond=None)[0]

    # Calculate genetic variance (Vg) and residual variance (Ve)
    Vg = np.var(K @ y_adjusted)
    Ve = np.var(y_adjusted - K @ y_adjusted)

    # Calculate heritability (h2)
    heritability = Vg / (Vg + Ve)
    heritability_percent = heritability * 100

    print(f"Heritability: {heritability:.4f}")
    print(f"Heritability (%): {heritability_percent:.2f}%")

    from sklearn.linear_model import LinearRegression

    # Fit linear model to calculate log likelihood
    lm = LinearRegression().fit(X_with_pcs, y)
    y_pred = lm.predict(X_with_pcs)
    residuals = y - y_pred

    # Calculate the residual sum of squares
    RSS = np.sum(residuals**2)

    # Calculate the log likelihood
    n = len(y)
    sigma2 = RSS / n
    log_likelihood = -0.5 * n * (np.log(2 * np.pi * sigma2) + 1)

    # Calculate -2 * log likelihood
    neg2_log_likelihood = -2 * log_likelihood

    print(f"-2 * Log Likelihood: {neg2_log_likelihood:.2f}")

    # Plot heritability
    plt.figure(figsize=(4, 4))
    plt.bar(['Heritability'], [heritability_percent], color='#8cc5e3')
    plt.ylabel('Heritability (%)')
    plt.title('Estimated Heritability')
    plt.suptitle(
        f'Heritability: {heritability:.4f}\nHeritability: {heritability:.2f}%\n-2 * Log Likelihood: {neg2_log_likelihood:.2f}',
        fontsize=12, y=1.05, ha='center')
    plt.ylim(0, 100)
    plt.tight_layout(rect=[0, 0, 1, 0.9])  # Adjust layout to make space for the title
    plt.show()
#h1()
h2()