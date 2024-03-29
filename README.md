# Trans-ancestry proteome-wide Mendelian randomization
In this study, we systematically estimated the causal role of 1,311 and 1,310 proteins, measured in populations from African and European ancestry respectively, on eight common diseases using a comprehensive ancestry-specific MR pipeline.
The results highlight the value of proteome-wide MR in informing the generalisability of drugs and drug targets across ancestries and illustrate the value of multi-cohort and biobank meta-analysis of genetic data for drug development.


We report our MR results in an openly accessible database: EpiGraphDB (https://epigraphdb.org/trans-ancestry-PWMR/). 

We also upload the discovery proteome-wide MR code here in `"R code"` folder and attach the example files in `"example file"` folder.

To start using the code, you need to install `TwoSampleMR`, `MendelianRandomization` and `ieugwasr` package.

```key
install.packages("remotes")
remotes::install_github("MRCIEU/TwoSampleMR")
```

```key
devtools::install_github("mrcieu/ieugwasr")
```

we include the following functions in R code for reference:
* LD clumping and F statistics in `instrument.R`
* steiger filtering, heterogenity and pleiotropy test in `2SMR.R`
* LDcheck in `LDcheck and colocalization.R`

The colocalization method: PWCOCO (https://github.com/jwr-git/pwcoco/). 

# Publication
Any other related information could be referenced from our paper: Proteome-wide Mendelian randomization in global biobank meta-analysis reveals multi-ancestry drug targets for common diseases (https://www.medrxiv.org/content/10.1101/2022.01.09.21268473v1)

# Data Sources
* GWAS data: [Global Biobank Meta-analysis Initiative (GBMI)](https://www.globalbiobankmeta.org/).
* pQTL data: [Zhang et al., 2021 BioRxiv](https://www.biorxiv.org/content/10.1101/2021.03.15.435533v1.full).


