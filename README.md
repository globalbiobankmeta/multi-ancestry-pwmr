# Effect of GLP1R siganlling and self-harm using Mendelian randomization
In this study, we estimated the effect of GLP1R signalling on self-harm outcomes using Mendelian randomization (MR), colocalization and polygenic score association. 
The results highlight that the weight loss effect of GLP1R signalling was likely to increase the risk of self-harm behaviours and raised concern about mental health safety of using GLP1R agonists for weight control, potentially with higher risk in women. 

We upload the discovery MR code here in `"Mendelian randomization"`.

To start using the code, you need to install `TwoSampleMR` and `ieugwasr` package.

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


