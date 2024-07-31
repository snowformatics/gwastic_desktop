import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from pysnptools.snpreader import Bed
from pysnptools.standardizer import Unit
from matplotlib.backends.backend_pdf import PdfPages

# File paths
bed_fn = "//filer-5/agruppen/BIT/lueck/gwas_test_data/gwastic_paper/small_set/example.bed"

# Read BED file
bed = Bed(bed_fn, count_A1=True).read()
print(f"Bed dimensions: {bed.iid_count} x {bed.sid_count}")

# Standardize the bed data (mean=0, sd=1)
bed_standardized = bed.standardize(Unit())

# Perform PCA
pca = PCA(n_components=10)
pca_result = pca.fit_transform(bed_standardized.val)

# Explained variance ratio
explained_variance_ratio = pca.explained_variance_ratio_

# Create a PDF file to save the plots
with PdfPages('PCA_plots.pdf') as pdf:
    # Scatter plots for the first four principal components
    for i in range(4):
        for j in range(i + 1, 4):
            plt.figure(figsize=(10, 8))
            plt.scatter(pca_result[:, i], pca_result[:, j], edgecolor='k', s=50)
            plt.xlabel(f'PC{i + 1}')
            plt.ylabel(f'PC{j + 1}')
            plt.title(f'PC{i + 1} vs PC{j + 1}')
            pdf.savefig()
            plt.close()

    # Scree plot (Explained Variance) for the first 10 components
    plt.figure(figsize=(10, 8))
    components = np.arange(1, 11)
    plt.bar(components, explained_variance_ratio * 100, color='grey')
    plt.xlabel('Principal Components')
    plt.ylabel('Variance Explained (%)')
    plt.title('Scree Plot')
    plt.xticks(components)
    pdf.savefig()
    plt.close()

print("PCA plots saved to PCA_plots.pdf")
