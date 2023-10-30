# About GWAStic

GWAStic is a software for Genome-Wide Association Study (GWAS) that combines traditional statistical methods with the power of Artificial Intelligence (AI) for comprehensive genetic analysis.

![ALT TEXT](https://github.com/snowformatics/gwastic_desktop/blob/e0743e1f67e5300d083a96441bbf505b5d7a7696/gwastic_desktop/images/gui.PNG)

## Table of Contents  
- [Installation](#1-installation)  
- [References](#2-references)  
- [Documentation](#3-documentation)  


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

We recommend to install Anaconda:
1. Install Anaconda (https://www.anaconda.com/distribution/)
2. `conda create --name gwastic_env python=3.9`
3. `conda activate gwastic_env`
4. `conda install pip`

Then install GWAStic via pip:

5. `pip install gwastic_desktop`

6. Open the Anaconda prompt and activate your GWAStic environment.<br/>`conda activate gwastic_env`<br/>

7. Type `gwastic` in the command line to start the software.

### Linux:

We recommend to install Anaconda:
1. `wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh`
2. `chmod +x Miniconda3-latest-Linux-x86_64.sh`
3. `./Miniconda3-latest-Linux-x86_64.sh`
4. `cd /home/username/miniconda3`
5. `source ~/miniconda3/bin/activate`

Then install GWAStic via pip:

6. `pip install gwastic_desktop`

7. `apt-get install libgl1`

8. Type `gwastic` in the command line to start the software.


# Supported input file formats
- VCF file format (including vcf.gz) and Plink BED (binary) format are supported for all GWAS methods. In case of vcf, you first must convert the genotype data to bed file format. 

[VCF example file](https://github.com/snowformatics/data/blob/cd8ac371fe669711430a6a4d7c00960082b3cd4b/gwastic_test_data/example.vcf.gz)

- Phenotypic data must be three columns (Family ID; Within-family ID; Value) text or CSV file delimited by *space*.

[Phenotype example file](https://github.com/snowformatics/data/blob/cd8ac371fe669711430a6a4d7c00960082b3cd4b/gwastic_test_data/pheno.csv)

# 2. References

# 3. Documentation





