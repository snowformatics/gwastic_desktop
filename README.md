# About GWAStic

GWAStic is a software for Genome-Wide Association Study (GWAS) that combines traditional statistical methods with the power of Artificial Intelligence (AI) for comprehensive genetic analysis.

![ALT TEXT](https://github.com/snowformatics/gwastic_desktop/blob/e0743e1f67e5300d083a96441bbf505b5d7a7696/gwastic_desktop/images/gui.PNG)

## Table of Contents  
- [Installation](#1-installation) 
- [Documentation](https://snowformatics.gitbook.io/product-docs/)  
- [References](#2-references)  
- [Acknowledgment](#3-acknowledgment)  


Key Features:
- Cross Platform 

- Comprehensive Genetic Analysis: GWAStic offers a wide range of methods to analyze your genomic data, allowing you to explore the associations between genetic variants and traits of interest comprehensively.

- AI-Enhanced Data Analysis: Harness the capabilities of machine learning and AI to uncover subtle patterns, interactions, and associations that may be missed by conventional statistical methods. 

- Genomic Prediction: Take your research to the next level by using GWAStic's advanced AI models for genomic prediction. Predict future health outcomes, disease risks, or phenotypic traits based on your genetic data and environmental factors.
- User-Friendly Interface: GWAStic's intuitive interface makes it accessible to both novice and experienced researchers. Seamlessly navigate through your data, perform analyses, and visualize results with ease.

- Customizable Workflows: Tailor your analysis to your specific research goals with customizable workflows. Define your parameters, select the appropriate statistical models, and integrate AI components as needed for a personalized analysis experience.

- Collaborative Research: Collaborate seamlessly with colleagues and share your findings securely within the platform. 

- Frequent Updates: Stay at the forefront of genetic research with regular software updates. GWAStic incorporates the latest advancements in GWAS and AI methodologies to keep your analyses up-to-date.

![myfile](https://github.com/snowformatics/gwastic_desktop/blob/08383abc5a0ba7920a542257b058094b85fe4446/gwastic_desktop/images/gui.gif)



# 1. Installation  

GWAStic software was build and successfully tested on Windows operating system (Windows 7 and 10).

### Windows:

> [!TIP]
> We recommend to install Anaconda and  for managing dependencies, it is often recommended to create a new environment for your project:
Install Anaconda (https://www.anaconda.com/distribution/)

```
conda create --name gwastic_env python=3.9
conda activate gwastic_env
conda install pip
```

> [!IMPORTANT]
> Install GWAStic via pip:
> 
`pip install gwastic_desktop`

> [!IMPORTANT]
> Run GWAStic:

Type `gwastic` in the Anaconda command line to start the software.

### Linux:

> [!TIP]
> We recommend to install Anaconda and  for managing dependencies, it is often recommended to create a new environment for your project:
`wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh`
`chmod +x Miniconda3-latest-Linux-x86_64.sh`
`./Miniconda3-latest-Linux-x86_64.sh`
`cd /home/username/miniconda3`
`source ~/miniconda3/bin/activate`
`conda create --name gwastic_env python=3.9`
`conda activate gwastic_env`

> [!IMPORTANT]
> Install GWAStic via pip:

8. `pip install gwastic_desktop`
> [!IMPORTANT]
> Run GWAStic:
9. Type `gwastic` in the command line to start the software.

### MacOS:

> [!TIP]
> We recommend to install Anaconda and  for managing dependencies, it is often recommended to create a new environment for your project:
1. Install Anaconda (https://www.anaconda.com/distribution/)

2. Open the downloaded .pkg file to launch the installer and follow the on-screen instructions.
3. Open Terminal. You can do this by pressing Cmd + Space to open Spotlight Search, typing "Terminal", and pressing Enter.
4. Create a new environment :
6. `conda create --name gwastic_env python=3.9`
7. `conda activate gwastic_env`

> [!IMPORTANT]
> Install GWAStic via pip:

`pip install gwastic_desktop`
9. For MacOS it's important to update matplotlib via conda:

`conda install matplotlib`


> [!IMPORTANT]
> Run GWAStic:
9. Type `gwastic` in the command line to start the software.



# Supported input file formats
- VCF file format (including vcf.gz) and Plink BED (binary) format are supported for all GWAS methods. In case of vcf, you first must convert the genotype data to bed file format. 

[VCF example file](https://github.com/snowformatics/data/blob/cd8ac371fe669711430a6a4d7c00960082b3cd4b/gwastic_test_data/example.vcf.gz)

- Phenotypic data must be three columns (Family ID; Within-family ID; Value) text or CSV file delimited by *space*.

[Phenotype example file](https://github.com/snowformatics/data/blob/cd8ac371fe669711430a6a4d7c00960082b3cd4b/gwastic_test_data/pheno.csv)

# 2. References

 2.1 - Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR, Bender D, Maller J, de Bakker PIW:
 Daly MJ & Sham PC (in press) PLINK: a toolset for whole-genome association and population-based linkage analysis. American Journal of Human Genetics.

 2.2 -  Lippert, C., Listgarten, J., Liu, Y. et al. FaST linear mixed models for genome-wide association studies. Nat Methods 8, 833–835 (2011). https://doi.org/10.1038/nmeth.1681

# 3. Acknowledgment
Gwastic has incorporated the FaST-LMM library (fastlmm.github.io), to enhance its Linear Mixed Models (LMM) feature. 
We thank Carl Kadie and David Heckerman for not only creating this exceptional tool but also providing outstanding support and discussions.




