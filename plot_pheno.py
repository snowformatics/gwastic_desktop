import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
#import imageio
import os

def pdf2():
    # Load the phenotype data from a space-delimited file
    file_path = 'pheno_gwas2.csv'  # Change this to your file path
    df = pd.read_csv(file_path, delim_whitespace=True, header=None, names=['ID1', 'ID2', 'Observation'])

    # Set a pastel color palette
    sns.set_palette("pastel")

    image_files = []

    # Create a PDF to save the plots and save individual images
    with PdfPages('test.pdf') as pdf:
        # Scatter plot: Observation vs individual
        plt.figure(figsize=(10, 6))
        plt.suptitle('Phenotype Data Distribution', fontsize=16)
        plt.scatter(df.index, df['Observation'], color='#8cc5e3')
        plt.xlabel('Individual')
        plt.ylabel('Observation')
        plt.title('Observation vs Individual')
        pdf.savefig()
        image_path = 'scatter_plot.png'
        plt.savefig(image_path)
        image_files.append(image_path)
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
        image_path = 'boxplot_violin_plot.png'
        plt.savefig(image_path)
        image_files.append(image_path)
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
        sns.kdeplot(df['Observation'], shade=True, color='#b8b8b8')
        plt.xlabel('Observation')
        plt.ylabel('Density')
        plt.title('Density Plot: Density vs Observation')

        # Save the combined plot to the PDF
        pdf.savefig()
        image_path = 'histogram_density_plot.png'
        plt.savefig(image_path)
        image_files.append(image_path)
        plt.close()

        # Cumulative density plot: Accumulative Density vs Observation
        plt.figure(figsize=(10, 6))
        sns.ecdfplot(df['Observation'], color='#8cc5e3')
        plt.xlabel('Observation')
        plt.ylabel('Cumulative Density')
        plt.title('Cumulative Density Plot: Accumulative Density vs Observation')
        pdf.savefig()
        image_path = 'cumulative_density_plot.png'
        plt.savefig(image_path)
        image_files.append(image_path)
        plt.close()

    # # Create GIF from the saved images
    # images = [imageio.imread(image_file) for image_file in image_files]
    # imageio.mimsave('plots.gif', images, duration=1)  # duration is in seconds
    #
    # # Clean up image files
    # for image_file in image_files:
    #     os.remove(image_file)

pdf2()



# import pandas as pd
# import matplotlib.pyplot as plt
# import seaborn as sns
# from matplotlib.backends.backend_pdf import PdfPages
#
# def pdf2():
#     import pandas as pd
#     import matplotlib.pyplot as plt
#     import seaborn as sns
#     from matplotlib.backends.backend_pdf import PdfPages
#
#     # Load the phenotype data from a space-delimited file
#     file_path = 'pheno_gwas2.csv'  # Change this to your file path
#     df = pd.read_csv(file_path, delim_whitespace=True, header=None, names=['ID1', 'ID2', 'Observation'])
#
#     # Set a grey color palette
#     # Set a pastel color palette
#     sns.set_palette("pastel")
#
#     # Create a PDF to save the plots
#     with PdfPages('test.pdf') as pdf:
#         # Scatter plot: Observation vs individual
#         plt.figure(figsize=(10, 6))
#         plt.suptitle('Phenotype Data Distribution', fontsize=16)
#         plt.scatter(df.index, df['Observation'], color='#8cc5e3')
#         plt.xlabel('Individual')
#         plt.ylabel('Observation')
#         plt.title('Observation vs Individual')
#         pdf.savefig()
#         plt.close()
#
#         # Boxplot and Violin plot side by side
#         plt.figure(figsize=(10, 6))
#
#         # Boxplot
#         plt.subplot(1, 2, 1)
#         sns.boxplot(y=df['Observation'], color='#8cc5e3')
#         plt.xlabel('Category')  # Add x-axis label
#         plt.title('Boxplot of Observations')
#
#         # Violin plot
#         plt.subplot(1, 2, 2)
#         sns.violinplot(y=df['Observation'], color='#b8b8b8')
#         plt.xlabel('Category')  # Add x-axis label
#         plt.title('Violin Plot of Observations')
#
#         # Save the combined plot to the PDF
#         pdf.savefig()
#         plt.close()
#
#         # Histogram and Density plot side by side
#         plt.figure(figsize=(12, 6))
#
#         # Histogram: Frequency vs Observation
#         plt.subplot(1, 2, 1)
#         plt.hist(df['Observation'], bins=30, color='#8cc5e3')
#         plt.xlabel('Observation')
#         plt.ylabel('Frequency')
#         plt.title('Histogram: Frequency vs Observation')
#
#         # Density plot: Density vs Observation
#         plt.subplot(1, 2, 2)
#         sns.kdeplot(df['Observation'], shade=True, color='#b8b8b8')
#         plt.xlabel('Observation')
#         plt.ylabel('Density')
#         plt.title('Density Plot: Density vs Observation')
#
#         # Save the combined plot to the PDF
#         pdf.savefig()
#         plt.close()
#
#         # Cumulative density plot: Accumulative Density vs Observation
#         plt.figure(figsize=(10, 6))
#         sns.ecdfplot(df['Observation'], color='#8cc5e3')
#         plt.xlabel('Observation')
#         plt.ylabel('Cumulative Density')
#         plt.title('Cumulative Density Plot: Accumulative Density vs Observation')
#         pdf.savefig()
#         plt.close()
#
#
#
#
# pdf2()