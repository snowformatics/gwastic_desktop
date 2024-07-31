import numpy as np
from pysnptools.snpreader import Bed, Pheno
import pysnptools.util as pstutil
from fastlmm.association import single_snp
#from fastlmm.util.util import merge_some
from pysnptools.standardizer import Unit
from fastlmm.inference import FastLMM
from pysnptools.kernelreader import SnpKernel

# File paths
#bed_fn = "//filer-5/agruppen/BIT/lueck/gwas_test_data/gwastic_paper/small_set/example.bed"
#pheno_fn = "pheno_gwas.csv"
#bed_fn = "//filer-5/agruppen/BIT/lueck/gwas_test_data/gwastic_paper/WGS300_005_0020.bed"
#pheno_file = "//filer-5/agruppen/BIT/lueck/NH_Paper/barley/phenoytpes/random_blues.txt"


# Read BED file
bed = Bed(bed_fn, count_A1=True).read()
print(f"Bed dimensions: {bed.iid_count} x {bed.sid_count}")

# Read phenotype file
pheno = Pheno(pheno_fn).read()

# Intersect bed and phenotype data
bed, pheno = pstutil.intersect_apply([bed, pheno])
print(f"Bed dimensions: {bed.iid_count} x {bed.sid_count}")

# Standardize the bed data
bed_standardized = bed.standardize(Unit())
print(f"Sample bed values: \n{bed_standardized.val[0:5,0:5]}")

# Create kinship matrix
kernel = SnpKernel(bed_standardized, standardizer=Unit())
kernel = kernel.read().standardize()  # defaults to DiagKtoN standardize
K = kernel.read().val
print(f"Sample kernel values: \n{K[0:5,0:5]}")

# Ensure the phenotype values are correctly formatted
y = pheno.val.flatten()

# Run Fast-LMM single_snp to get heritability estimates
result = single_snp(test_snps=bed_standardized, pheno=pheno, K0=K, leave_out_one_chrom=False, output_file_name=None)

# Extract heritability estimates (Nullh2)
heritability_estimates = result['Nullh2']
mean_heritability = np.mean(heritability_estimates)
heritability_percent = mean_heritability * 100

print(f"Mean Heritability: {mean_heritability:.4f}")
print(f"Mean Heritability (%): {heritability_percent:.2f}%")

# Optional: Print detailed heritability results
print(result[['SNP', 'Nullh2']])