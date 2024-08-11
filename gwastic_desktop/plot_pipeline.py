import pandas as pd
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np

import numpy as np

from pysnptools.snpreader import Bed, Pheno
import pysnptools.util as pstutil
from pysnptools.kernelreader import SnpKernel
from pysnptools.standardizer import Unit
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.decomposition import PCA
import pandas as pd

plt.switch_backend('Agg')


class Plot:
    """Plotting class."""

    def __init__(self):
        pass

    def plot_pheno_statistics(self, pheno_file, pheno_stats_name):
        #print (pheno.iid, pheno.val)
        # Load the phenotype data from a space-delimited file
        file_path = pheno_file  # Change this to your file path

        df = pd.read_csv(file_path, delim_whitespace=True, header=None, names=['ID1', 'ID2', 'Observation'])

        # Create a PDF to save the plots
        with PdfPages(pheno_stats_name) as pdf:
            # Scatter plot: Observation vs individual
            plt.figure(figsize=(10, 6))
            plt.scatter(df.index, df['Observation'], color='#8cc5e3')
            plt.xlabel('Individual')
            plt.ylabel('Observation')
            plt.title('Observation vs Individual')
            pdf.savefig()
            plt.close()

            # Boxplot and Violin plot side by side
            plt.figure(figsize=(10, 6))

            # Boxplot
            plt.subplot(1, 2, 1)
            sns.boxplot(y=df['Observation'], color='#8cc5e3')
            plt.xlabel('Category')  # Add x-axis label
            plt.title('Boxplot of Observations')

            # Violin plot
            plt.subplot(1, 2, 2)
            sns.violinplot(y=df['Observation'], color='#b8b8b8')
            plt.xlabel('Category')  # Add x-axis label
            plt.title('Violin Plot of Observations')

            # Save the combined plot to the PDF
            pdf.savefig()
            plt.close()

            # Histogram and Density plot side by side
            plt.figure(figsize=(12, 6))

            # Histogram: Frequency vs Observation
            plt.subplot(1, 2, 1)
            plt.hist(df['Observation'], bins=30, color='#8cc5e3')
            plt.xlabel('Observation')
            plt.ylabel('Frequency')
            plt.title('Histogram: Frequency vs Observation')

            # Density plot: Density vs Observation
            plt.subplot(1, 2, 2)
            sns.kdeplot(df['Observation'], fill=True, color='#b8b8b8')
            plt.xlabel('Observation')
            plt.ylabel('Density')
            plt.title('Density Plot: Density vs Observation')

            # Save the combined plot to the PDF
            pdf.savefig()
            plt.close()

            # Cumulative density plot: Accumulative Density vs Observation
            plt.figure(figsize=(10, 6))
            sns.ecdfplot(df['Observation'], color='#8cc5e3')
            plt.xlabel('Observation')
            plt.ylabel('Cumulative Density')
            plt.title('Cumulative Density Plot: Accumulative Density vs Observation')
            pdf.savefig()
            plt.close()

    def plot_geno_statistics(self, bed, pheno, geno_stats_name):

        bed = bed.read()
        pheno = pheno.read()

        # Calculate MAF (Minor Allele Frequency)
        allele_sums = bed.val.sum(axis=0)  # Sum of allele counts for each SNP
        allele_counts = 2 * bed.iid_count  # Each individual has 2 alleles
        allele_frequencies = allele_sums / allele_counts
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

        # Kinship
        kernel2 = bed_standardized.val @ bed_standardized.val.T / bed_standardized.sid_count
        sample_names = [''.join(iid) for iid in bed.iid[:, 1]]

        # Hertitabilty
        # Add principal components as covariates
        X_with_pcs = np.hstack((pca_result, np.ones((pca_result.shape[0], 1))))  # Include an intercept

        # Adjust phenotype for population structure
        pheno = pheno.read()
        y = pheno.val.flatten()
        y_adjusted = y - X_with_pcs @ np.linalg.lstsq(X_with_pcs, y, rcond=None)[0]

        # Calculate genetic variance (Vg) and residual variance (Ve)
        Vg = np.var(kernel2 @ y_adjusted)
        Ve = np.var(y_adjusted - kernel2 @ y_adjusted)

        # Calculate heritability (h2)
        heritability = Vg / (Vg + Ve)
        heritability_percent = heritability * 100

        # Create a PDF to save the plots
        with PdfPages(geno_stats_name) as pdf:
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

            pdf.savefig()
            plt.close()

            plt.figure(figsize=(12, 4))
            plt.suptitle('Heritability and Explained Variance Ratio', fontsize=16)
            plt.subplots_adjust(top=0.75, wspace=0.5)

            # Heritability percentage barchart
            plt.subplot(1, 2, 1)
            plt.bar(['Heritability'], [heritability_percent], color='#8cc5e3')
            plt.ylim(0, 100)
            #plt.xlim(0, 2)
            plt.ylabel('Percentage')
            plt.title('Heritability (%) ' + str(heritability_percent))

            # Explained variance ratio plot
            plt.subplot(1, 2, 2)
            plt.bar([f'PC{i + 1}' for i in range(len(explained_variance_ratio))], explained_variance_ratio,
                    color='#b8b8b8')
            plt.ylabel('Explained Variance Ratio')
            plt.title('Explained Variance Ratio by Principal Components')

            # Save the combined plots to the PDF
            pdf.savefig()
            plt.close()

            try:
                kinship_df = pd.DataFrame(kernel2, index=sample_names, columns=sample_names)
                g = sns.clustermap(kinship_df, cmap='vlag', figsize=(15, 15), xticklabels=1, yticklabels=1,
                                   method='single', metric='euclidean')

                # Adjust the title positioning
                plt.subplots_adjust(top=0.95)
                g.fig.suptitle('Hierarchically-Clustered Kinship Matrix', fontsize=16)

                # Adjust the font size for labels
                plt.setp(g.ax_heatmap.get_xticklabels(), fontsize=4, rotation=90)
                plt.setp(g.ax_heatmap.get_yticklabels(), fontsize=4, rotation=0)
                x0, _y0, _w, _h = g.cbar_pos
                g.ax_cbar.set_position([x0, 0.9, g.ax_row_dendrogram.get_position().width - 0.1, 0.05])
                # Save the heatmap to the PDF
                pdf.savefig(g.fig)
                plt.close()
            except RecursionError:
                pass

