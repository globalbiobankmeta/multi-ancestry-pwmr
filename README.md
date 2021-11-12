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

Th colocalization method we choose is PWCOCO (https://github.com/jwr-git/pwcoco/). 

Any other related information could be referenced from our paper (https://outlook.office.com/mail/group/groups.bristol.ac.uk/grp-trans-ethnicpwas-gbmi/files/sxs/sp/701F7E1B-F49A-421F-855C-B72A1CC39EB4/). 
