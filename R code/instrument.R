library(TwoSampleMR)
library(ieugwasr)

eapqtl<- read.table(file="eapqtl.txt", header=T, sep="\t", quote="", stringsAsFactors=F, check.names=F)

exposure_samplesize <- 7213
outcome_samplesize <- 1219215
ncase <- 100736
ncontrol <- 1118479
prevalence <- ncase/outcome_samplesize

###format data
tmp <- data.frame(
  SNP = eapqtl$ID,
  beta = eapqtl$BETA,
  se = eapqtl$SE,
  effect_allele = eapqtl$ea,
  eaf = eapqtl$A1_FREQ,
  other_allele = eapqtl$oa,
  phenotype_id = eapqtl$`SOMAmer ID`,
  pval = eapqtl$P,
  chr = eapqtl$CHROM,
  pos = eapqtl$POS,
  samplesize = exposure_samplesize
)
exposure_dat <- format_data(tmp, type="exposure", effect_allele_col = "effect_allele",other_allele_col = "other_allele",phenotype_col = "phenotype_id",eaf_col = "eaf",samplesize_col = "samplesize",chr_col = "chr",pos_col = "pos")

#LD clumping
ld_reflookup <- function(rsid, pop='EUR')
{
  res <- api_query('ld/reflookup',
                   query = list(
                     rsid = rsid,
                     pop = pop
                   )
  ) %>% get_query_content()
  if(length(res) == 0)
  {
    res <- character(0)
  }
  return(res)
}
refsnp<-ld_reflookup(exposure_dat$SNP)
exp1<-exposure_dat[exposure_dat$SNP %in% refsnp,]
exp2<-exposure_dat[!exposure_dat$SNP %in% refsnp,]
exp1 <-clump_data(exp1,clump_kb=10000,clump_r2=0.6,clump_p1=1,clump_p2=1,pop="EUR")
exposure_dat<-rbind(exp1,exp2)


###F-statistics for each SNPs
exposure_dat <- cbind(exposure_dat,fstatistics=1)
for (s in 1:nrow(exposure_dat)){
  z <- exposure_dat[s,2]/exposure_dat[s,3]
  pve <- z^2/(exposure_samplesize+z^2)
  exposure_dat[s,"fstatistics"] <- (exposure_samplesize-2)*pve/(1-pve)
}
exposure_dat <- exposure_dat[exposure_dat$fstatistics>10,]


