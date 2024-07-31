import numpy as np
from fastlmm.util import example_file # Download and return local file name
from pysnptools.snpreader import Bed, Pheno
import pysnptools.util as pstutil
from pysnptools.kernelreader import SnpKernel
from pysnptools.standardizer import Unit
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.decomposition import PCA
import pandas as pd

bed_fn = "//filer-5/agruppen/BIT/lueck/gwas_test_data/gwastic_paper/small_set/example.bed"
pheno_fn = "pheno_gwas.csv"
bed_fn = "//filer-5/agruppen/BIT/lueck/gwas_test_data/gwastic_paper/WGS300_005_0020.bed"
pheno_fn = "pheno_gwas2.csv"


bed = Bed(bed_fn, count_A1=True).read()


# Calculate allele frequencies
allele_sums = bed.val.sum(axis=0)  # Sum of allele counts for each SNP
allele_counts = 2 * bed.iid_count  # Each individual has 2 alleles
allele_frequencies = allele_sums / allele_counts
# Calculate MAF (Minor Allele Frequency)
maf = np.minimum(allele_frequencies, 1 - allele_frequencies)
# Calculate heterozygosity for each individual (iid)
heterozygosity = np.mean(bed.val == 1, axis=1)
# Calculate heterozygosity for each marker (sid)
marker_heterozygosity = np.mean(bed.val == 1, axis=0)

# Perform PCA
bed_standardized = bed.standardize(Unit())
pca = PCA(n_components=10)
pca_result = pca.fit_transform(bed_standardized.val)
# Explained variance ratio
explained_variance_ratio = pca.explained_variance_ratio_

# create kinship matrix
sample_names = [''.join(iid) for iid in bed.iid[:, 1]]
#kernel = SnpKernel(bed_standardized, standardizer=Unit())
#kernel = kernel.read().standardize()  # defaults to DiagKtoN standardize
#print(f"sample kernel values: \n{kernel.val[0:5,0:5]}")

# Behind the scenes this is the same as:
#kernel2 = snps * snps.T / snps.shape[1]
kernel2 = bed_standardized.val @ bed_standardized.val.T / bed_standardized.sid_count
#print(f"sample kernel2 values: \n{kernel2[0:5,0:5]}")

# Set a pastel color palette
#sns.set_palette("pastel")

# Create a PDF to save the plots
with PdfPages('combined_plots.pdf') as pdf:
    # Create a new page with statistics as text
    plt.figure(figsize=(18, 6))  # Standard letter size
    plt.axis('off')
    plt.text(0.5, 0.9, 'GWAS Summary Statistics', fontsize=16, ha='center')
    plt.text(0.1, 0.8, f'Number of SNPs: {bed.sid_count}', fontsize=12)
    plt.text(0.1, 0.75, f'Number of Samples: {bed.iid_count}', fontsize=12)
    plt.text(0.1, 0.7, f'Heritability: {0.72:.4f}', fontsize=12)
    plt.text(0.1, 0.65, f'Heritability (%): {72:.2f}%', fontsize=12)
    plt.text(0.1, 0.6, f'-2 * Log Likelihood: {400:.2f}', fontsize=12)
    plt.text(0.1, 0.55, f'Principal Components Explained Variance:', fontsize=12)
    for i, var_ratio in enumerate(explained_variance_ratio, start=1):
        plt.text(0.15, 0.55 - i * 0.03, f'PC{i}: {var_ratio:.4f}', fontsize=12)
    plt.tight_layout()
    pdf.savefig()
    plt.close()

    plt.figure(figsize=(18, 6))  # Adjust the figure size as needed
    plt.suptitle('Marker density\n', fontsize=16)

    # Histogram of Minor Allele Frequency (MAF)
    plt.subplot(1, 3, 1)
    sns.histplot(maf, bins=9, kde=False, color='#8cc5e3', edgecolor=None)
    plt.xlabel('Minor Allele Frequency (MAF)')
    plt.ylabel('Frequency')
    plt.title('Histogram of Minor Allele Frequency (MAF)')

    # Histogram of Individual Heterozygosity
    plt.subplot(1, 3, 2)
    sns.histplot(heterozygosity, bins=5, kde=False, color='#b8b8b8', edgecolor=None)
    plt.xlabel('Individual Heterozygosity')
    plt.ylabel('Frequency')
    plt.title('Histogram of Individual Heterozygosity')

    # Histogram of Marker Heterozygosity
    plt.subplot(1, 3, 3)
    sns.histplot(marker_heterozygosity, bins=10, kde=False, color='#8cc5e3', edgecolor=None)
    plt.xlabel('Marker Heterozygosity')
    plt.ylabel('Frequency')
    plt.title('Histogram of Marker Heterozygosity')

    # Save the combined plots to the PDF
    pdf.savefig()
    plt.savefig('Marker.png')
    plt.close()

    # Create a figure for PCA scatter plots
    plt.figure(figsize=(18, 12))  # Adjusted figure size to fit all plots on one page
    plt.suptitle('Principal component analysis', fontsize=16)
    # Define the hex color codes
    color1 = '#b8b8b8'
    color2 = '#8cc5e3'

    plot_number = 1
    for i in range(4):
        for j in range(i + 1, 4):
            plt.subplot(2, 3, plot_number)
            # Alternate colors based on the index
            color = color1 if plot_number % 2 == 0 else color2
            sns.scatterplot(x=pca_result[:, i], y=pca_result[:, j], edgecolor=None, s=50, color=color)
            plt.xlabel(f'PC{i + 1}')
            plt.ylabel(f'PC{j + 1}')
            plt.title(f'PC{i + 1} vs PC{j + 1}')
            plot_number += 1

    # Save the PCA scatter plots to the PDF
    plt.savefig('pca.png')
    pdf.savefig()
    plt.close()

    # Scree plot (Explained Variance) for the first 10 components on a new page
    # plt.figure(figsize=(10, 8))
    #
    # components = np.arange(1, 11)
    # plt.bar(components, explained_variance_ratio * 100, color=sns.color_palette("pastel")[3], edgecolor=None)
    #
    # plt.xlabel('Principal Components')
    # plt.ylabel('Variance Explained (%)')
    # plt.title('Scree Plot')
    # plt.xticks(components)
    #
    # # Save the scree plot to the PDF
    # pdf.savefig()
    # plt.close()


    # Create the hierarchically-clustered heatmap on a new page
    kinship_df = pd.DataFrame(kernel2, index=sample_names, columns=sample_names)
    g = sns.clustermap(kinship_df, cmap='vlag', figsize=(15, 15), xticklabels=1, yticklabels=1)

    # Adjust the title positioning
    plt.subplots_adjust(top=0.95)
    g.fig.suptitle('Hierarchically-Clustered Kinship Matrix', fontsize=16)

    # Adjust the font size for labels
    plt.setp(g.ax_heatmap.get_xticklabels(), fontsize=4, rotation=90)
    plt.setp(g.ax_heatmap.get_yticklabels(), fontsize=4, rotation=0)
    x0, _y0, _w, _h = g.cbar_pos
    g.ax_cbar.set_position([x0, 0.9, g.ax_row_dendrogram.get_position().width-0.1, 0.05])

    # Save the heatmap to the PDF
    pdf.savefig(g.fig)
    plt.savefig('kinship.png', dpi=300)
    plt.close()