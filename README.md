# DPW WES

This repository contains the code for performing analysis in the DPW_WES project.

## Table of Contents

- [System Requirements](#system-requirements)
- [Exome-Wide Association Study (ExWAS)](#exome-wide-association-study-exwas)
- [Burden Heritability Regression (BHR)](#burden-heritability-regression-bhr)
- [Expression Analysis](#expression-analysis)
- [License](#license)

## System Requirements

This package is supported for Linux. The package has been tested on the following systems:

Linux: x86_64 GNU/Linux

LSB Version:	:core-4.1-amd64:core-4.1-noarch

Distributor ID:	CentOS

Description:	CentOS Linux release 7.9.2009 (Core)

Release:	7.9.2009

Codename:	Core

R version: R/4.1.0

### R dependencies

The required R dependencies for each script are described within the respective R script files. Make sure to install the necessary packages before running the code.

## Exome-Wide Association Study (ExWAS)

In this project, we conducted an exome-wide association analysis using SAIGE GENE+. SAIGE GENE+ is a powerful tool for rare variant association tests. For detailed information about SAIGE GENE+, please refer to the publication by Zhou et al., "SAIGE-GENE+ improves the efficiency and accuracy of set-based rare variant association tests" ([Nature Genetics, 2022](https://doi.org/10.1038/s41588-022-01178-w)).

To get started with SAIGE GENE+, follow the installation instructions provided on the [GitHub repository](https://github.com/weizhouUMICH/SAIGE). Download the main R file from the repository.

For fitting null generalized linear mixed models (GLMMs) for binary traits (e.g., Alzheimer's disease case-control status), use the script `S1_binary.sh`. For fitting null GLMMs for quantitative traits (e.g., drinks per week), use `S1_quantitative.sh`.

To perform gene-based association tests, use `S2_Gene-based.sh`. For single-variant association tests, use `S2_singlevariant.sh`. Make sure to configure the parameters according to the [SAIGE documentation](https://saigegit.github.io/SAIGE-doc/).

To submit tasks to the cluster, you can use the provided example script `qsub_task.sh`. Adjust the script as needed and run it for each chromosome.

An example output of gene-based association tests is provided as `DPW_GeneBased.txt.gz`, and an example output of single-variant association tests is provided as `DPW_SingleRareAssoc.txt.gz`. The individual data used in this analysis were obtained from the [UK Biobank](https://www.ukbiobank.ac.uk/).

## Burden Heritability Regression (BHR)

We used Burden Heritability Regression (BHR) to estimate the heritability explained by ultra-rare and rare variants in each gene. For a detailed description of BHR, please refer to the publication by Weiner et al., "Polygenic architecture of rare coding variation across 394,783 exomes" ([Nature, 2023](https://doi.org/10.1038/s41586-022-05684-z)).

To install BHR, use the following command in R:

```R
devtools::install_github("ajaynadig/bhr")
```

Run the BHR.R script to initiate the BHR analysis:

```shell
$ Rscript BHR.R 
```

## Expression Analysis

### Tissue Enrichment Analysis

We performed tissue enrichment analysis using the TissueEnrich R package. For installation and usage instructions, please visit the [Bioconductor website](https://bioconductor.org/packages/release/bioc/html/TissueEnrich.html).

To start the tissue enrichment analysis, run the `TissueEnrich.R` script:

```shell
$ Rscript TissueEnrich.R
```

### Human Protein atlas data

The rna_tissue_hpa.tsv.zip file was downloaded from the Human Protein Atlas [website](https://www.proteinatlas.org/).

### Single cell RNA expression
To analyze single-cell RNA expression data, run the `scRNAliver.R` script:
```shell
$ Rscript scRNAliver.R 
```
## License
This project is covered under the MIT License.

