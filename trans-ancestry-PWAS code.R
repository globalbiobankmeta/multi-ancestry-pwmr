########discovery PWAS main analysis for Asthma####################################################
######European ancestry#########
rm(list=ls())
gc()

library(readxl)
library(TwoSampleMR)
library(MendelianRandomization)
library(ieugwasr)

exposure_samplesize <- 7213
outcome_samplesize <- 1219215
ncase <- 100736
ncontrol <- 1118479
prevalence <- ncase/outcome_samplesize

ea <- read_excel("/Users/tq20202/Desktop/project2/asthma/ea/conditional_ea.xlsx")

fname<-list.files("/Users/tq20202/Desktop/project2/ea_summary")
id=strsplit(fname[1],"[.]")[[1]][2]
name<-cbind(filename=fname[1],id)
for ( i in 2: length(fname)){
  temp<-cbind(fname[i],strsplit(fname[i],"[.]")[[1]][2])
  name<-rbind(name,temp)
}
name<-as.data.frame(name)

setwd("/Users/tq20202/Desktop/project2/ea_summary")
eares<-c()
for (i in 1: nrow(ea)){
  f<-name[name$id==as.character(ea[i,1]),1]
  if(length(f)==0){pi=rep(as.character(ea[i,1]),14)}
  else{temp <- read.table(file=f, header=T, sep="\t", quote="", stringsAsFactors=F, check.names=F)
      temp<-rbind(colnames(temp),temp[1:nrow(temp),])
      colnames(temp)<-c("CHROM","POS","ID",	"REF","ALT","A1",	"A1_FREQ","TEST",	"OBS_CT",	"BETA",	"SE",	"T_STAT",	"P",	"ERRCODE")
      pi<-temp[temp$ID==ea[i,10]$`Conditional independent cis-SNP`,]}
  eares<-rbind(eares,pi)
}
write.table(eares, file="/Users/tq20202/Desktop/project2/asthma/ea/eares.txt", row.names=F, col.names=T, sep="\t", quote=F)

eares<-cbind(ea[,1:5],eares)
eapqtl<-eares[eares$TEST=="ADD",] ##match summary stats
eapqtl<-cbind(eapqtl[,1:12],factor=0,ea=0,oa=0,eapqtl[,13:18])
eapqtl$A1_FREQ<-as.numeric(eapqtl$A1_FREQ)
for (i in 1:nrow(eapqtl)){
  ref<-eapqtl[i,"REF"]
  alt<-eapqtl[i,"ALT"]
  if (eapqtl[i,"A1"]==ref){eapqtl[i,"factor"]=-1;eapqtl[i,"ea"]=eapqtl[i,"A1"];eapqtl[i,"oa"]=eapqtl[i,"ALT"]}
  if (eapqtl[i,"A1"]==alt){eapqtl[i,"factor"]=1;eapqtl[i,"ea"]=eapqtl[i,"A1"];eapqtl[i,"oa"]=eapqtl[i,"REF"]}
}
unique(eapqtl$factor)

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
write.table(exp1, file="/Users/tq20202/Desktop/project2/asthma/ea/exp1.txt", row.names=F, col.names=T, sep="\t", quote=F)
write.table(exp2, file="/Users/tq20202/Desktop/project2/asthma/ea/exp2.txt", row.names=F, col.names=T, sep="\t", quote=F)

###F-statistics for each SNPs
exposure_dat <- cbind(exposure_dat,fstatistics=1)
for (s in 1:nrow(exposure_dat)){
  z <- exposure_dat[s,2]/exposure_dat[s,3]
  pve <- z^2/(exposure_samplesize+z^2)
  exposure_dat[s,"fstatistics"] <- (exposure_samplesize-2)*pve/(1-pve)
}
print(min(exposure_dat$fstatistics))
print(max(exposure_dat$fstatistics))
write.table(exposure_dat, file="/Users/tq20202/Desktop/project2/asthma/ea/exposure_dat_ea.txt", row.names=F, col.names=T, sep="\t", quote=F)

exposure_dat <- exposure_dat[exposure_dat$fstatistics>10,]
length(unique(exposure_dat$SNP))
length(unique(exposure_dat$exposure))

##outcome
Asthma_Bothsex_nfe_inv_var_meta_GBMI_052021 <- read.delim("/Users/tq20202/outcomedata/Asthma_Bothsex_nfe_inv_var_meta_GBMI_052021.txt", header=FALSE, comment.char="#")
outdate <- Asthma_Bothsex_nfe_inv_var_meta_GBMI_052021[Asthma_Bothsex_nfe_inv_var_meta_GBMI_052021$V5 %in% exposure_dat$SNP,]
length(unique(outdate$V5))
nosnp<-setdiff(exposure_dat$SNP,outdate$V5)
mexp<-exposure_dat[exposure_dat$SNP %in% nosnp,c("SNP","effect_allele.exposure","eaf.exposure","other_allele.exposure","beta.exposure")]
mexp<-mexp[!duplicated(mexp$SNP),]
######obtain unmatched snp outcome with proxy snp
eaproxy <- read.table(file="/Users/tq20202/Desktop/project2/EUR_pQTL_proxy.txt", header=T, sep="\t", quote="", stringsAsFactors=F, check.names=F)
eaproxy<-eaproxy[eaproxy$snp!=eaproxy$proxy,]
eaproxy<-eaproxy[eaproxy$snp %in% nosnp,]
s<-unique(eaproxy$snp)
proxydata<-c()
for (i in 1:length(s)){
  snp=s[i]
  d<-c()
  d2<-c()
  me<-mexp[mexp$SNP==snp,"effect_allele.exposure"]
  mef<-mexp[mexp$SNP==snp,"eaf.exposure"]
  mo<-mexp[mexp$SNP==snp,"other_allele.exposure"]
  temp<-eaproxy[eaproxy$snp==snp,"proxy"]
  out<-Asthma_Bothsex_nfe_inv_var_meta_GBMI_052021[Asthma_Bothsex_nfe_inv_var_meta_GBMI_052021$V5 %in% temp,]
  if (nrow(out)==0){
    d<-rep(NA,18)
  } else if (me %in% out$V4 & min(abs(out[out$V4==me,"V6"]-mef))<0.1){
    temp<-out[out$V4==me,]
    minf=min(abs(temp$V6-mef))
    d<-temp[abs(temp$V6-mef)==minf,]
    d<-cbind(d,snp)
  } else if (mo %in% out$V4 & min(abs((1-out[out$V4==mo,"V6"])-mef))<0.1){
    temp<-out[out$V4==mo,]
    minf=min(abs((1-temp$V6)-mef))
    d2<-temp[abs((1-temp$V6)-mef)==minf,]
    d2$V7=-d2$V7
    d2<-cbind(d2,snp)
  } else {
    minf=min(abs(out$V6-mef))
    d<-out[abs(out$V6-mef)==minf,]
    d<-cbind(d,snp)
  }
proxydata<-rbind(proxydata,d,d2)
}
proxy<-na.omit(proxydata[,c(1:9,18)])
msnp<-c()
for (i in 1:nrow(proxy)){
  temp=proxy[i,10]
  temp2<-mexp[mexp$SNP==temp,2:5]
  msnp<-rbind(msnp,temp2)
}
proxy<-cbind(proxy,msnp)
proxy<-proxy[abs(proxy$V6-proxy$eaf.exposure)<0.1,]
length(unique(proxy$snp))
write.table(proxy, file="/Users/tq20202/Desktop/project2/asthma/ea/proxy_ea.txt", row.names=F, col.names=T, sep="\t", quote=F)

proxy<-proxy[!duplicated(proxy$snp),]
proxy_data<-data.frame(SNP=proxy$snp,
                       effect_allele=proxy$effect_allele.exposure,
                       eaf=proxy$eaf.exposure,
                       other_allele=proxy$other_allele.exposure,
                       beta=proxy$V7,
                       se=proxy$V8,
                       pval=proxy$V9)
out_data<-data.frame(SNP=outdate$V5,
                     effect_allele=outdate$V4,
                     eaf=outdate$V6,
                     other_allele=outdate$V3,
                     beta=outdate$V7,
                     se=outdate$V8,
                     pval=outdate$V9)
outcome_dat <- rbind(out_data,proxy_data)
length(unique(outcome_dat$SNP))
outcome_dat  <- format_data(outcome_dat , type="outcome")
write.table(outcome_dat, file="/Users/tq20202/Desktop/project2/asthma/ea/outcome_dat_ea.txt", row.names=F, col.names=T, sep="\t", quote=F)

#####2SMR for matched data#################################################
inter<-intersect(exposure_dat$SNP,outcome_dat$SNP)
exposure_dat <- exposure_dat[exposure_dat$SNP %in% inter,]
outcome_dat <- outcome_dat[outcome_dat$SNP %in% inter,]
harmdat <- harmonise_data(exposure_dat, outcome_dat, action=2)
sum(abs(harmdat$eaf.exposure-harmdat$eaf.outcome)<0.1)
tharm=harmdat
tharm$test = tharm$effect_allele.exposure == tharm$effect_allele.outcome
nrow(tharm[tharm$test==FALSE,])
harmdat$mr_keep<-TRUE  #consider ambiguous palindromic SNPs (for the "MendelianRandomization" package)

test<-harmdat
test$new<- paste(test$SNP,test$exposure,sep="_")
test<-test[,-24]
colnames(test)[31]<-"exposure"
test<-test[,-26]
test <- format_data(test, type="exposure", snp_col = "SNP", pval_col = "pval.exposure",
                    beta_col = "beta.exposure", se_col = "se.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure",
                    phenotype_col = "exposure",eaf_col = "eaf.exposure",    
                    samplesize_col = "samplesize.exposure")
testout <- rbind(out_data,proxy_data)
testout  <- format_data(testout , type="outcome")
inter<-intersect(test$SNP,testout$SNP)
test <- test[test$SNP %in% inter,]
testout <- testout[testout$SNP %in% inter,]
testharm <- harmonise_data(test, testout, action=2)
testharm$mr_keep<-TRUE
testres<-mr(testharm)
###steiger filtering
testres <- generate_odds_ratios(testres)
testres <- cbind(SNP=1,testres,beta.exp=1,beta.out=1,se.exp=1,se.out=1,eaf.exp=1,eaf.out=1,p.exp=1,p.out=1)
for (r in 1: nrow(testharm)){
  testres[r,1] <- testharm[testharm$id.exposure==testres[r,"id.exposure"],"SNP"]
  testres[r,16:23] <- testharm[testharm$id.exposure==testres[r,"id.exposure"],c("beta.exposure","beta.outcome","se.exposure",
                                                                                "se.outcome","eaf.exposure","eaf.outcome",
                                                                                "pval.exposure","pval.outcome")]
}
testres <- na.omit(testres)
testres <- cbind(testres,rsq.exposure=1,rsq.outcome=1,samplesize.exposure=exposure_samplesize,samplesize.outcome=outcome_samplesize)
testres$rsq.exposure <- (get_r_from_pn(p=testres$p.exp, n=exposure_samplesize))^2
testres$rsq.outcome <- (get_r_from_lor(lor=log(testres$or),af=testres$eaf.out,ncase=ncase,ncontrol=ncontrol,prevalence=prevalence))^2
st <- psych::r.test(
  n = testres$samplesize.exposure, 
  n2 = testres$samplesize.outcome, 
  r12 = sqrt(testres$rsq.exposure), 
  r34 = sqrt(testres$rsq.outcome)
)
testres$steiger_dir <- testres$rsq.exposure > testres$rsq.outcome
testres$steiger_pval <- pnorm(-abs(st$z)) * 2

harmdat$match<-paste(harmdat$SNP,harmdat$exposure,sep="_")
harmdat$steiger_dir <- 0
harmdat$steiger_pval <- 0
for (i in 1:nrow(harmdat)){
  harmdat[i,33] <- testres[testres$exposure==harmdat[i,"match"],"steiger_dir"]
  harmdat[i,34] <- testres[testres$exposure==harmdat[i,"match"],"steiger_pval"]
}
harmdat<-harmdat[,-32]
harmdat<-harmdat[harmdat$steiger_dir==1,]
length(unique(harmdat$exposure))
write.table(harmdat, file="/Users/tq20202/Desktop/project2/asthma/ea/harmdat_ea.txt", row.names=F, col.names=T, sep="\t", quote=F)

#Generalized IVW  and MR-Egger for the "MendelianRandomization" package
exp<-unique(harmdat$exposure)
mrres <- c()
for (i in 1:length(exp)){
  dat <- harmdat[harmdat$exposure==exp[i],]
  if (nrow(dat)==1){
    result <- mr_wald_ratio(b_exp= dat$beta.exposure, b_out=dat$beta.outcome, se_exp= dat$se.exposure, se_out= dat$se.outcome)
    result <- data.frame(id=exp[i],method="Wald_ratio",nsnp=result$nsnp,b=result$b,se=result$se,CIlower=NA,CIupper=NA,
                         pval=result$pval,intercept=NA,intercept_se=NA,inter_CIlower=NA,inter_CIupper=NA,
                         intercept_pval=NA,hetero_Q=NA,hetero_pvale=NA)}
  else if (nrow(dat)==2 & all(dat$SNP %in% exp1$SNP)){rho <- ld_matrix(dat$SNP, pop = "EUR")
  result <- MendelianRandomization::mr_ivw(mr_input(bx = dat$beta.exposure,
                                                    bxse = dat$se.exposure,
                                                    by = dat$beta.outcome,
                                                    byse = dat$se.outcome,
                                                    cor = rho))
  result <- data.frame(id=exp[i],
                       method="IVW",
                       nsnp=result$SNPs,
                       b=result$Estimate,
                       se=result$StdError,
                       CIlower=result$CILower,
                       CIupper=result$CIUpper,
                       pval=result$Pvalue,
                       intercept=NA,
                       intercept_se=NA,
                       inter_CIlower=NA,
                       inter_CIupper=NA,
                       intercept_pval=NA,
                       hetero_Q=result$Heter.Stat[1],
                       hetero_pvale=result$Heter.Stat[2])}
  else if (nrow(dat)==2 & any(dat$SNP %in% exp2$SNP)){result <- MendelianRandomization::mr_ivw(mr_input(bx = dat$beta.exposure,
                                                                                                        bxse = dat$se.exposure,
                                                                                                        by = dat$beta.outcome,
                                                                                                        byse = dat$se.outcome))
  result <- data.frame(id=exp[i],
                       method="IVW",
                       nsnp=result$SNPs,
                       b=result$Estimate,
                       se=result$StdError,
                       CIlower=result$CILower,
                       CIupper=result$CIUpper,
                       pval=result$Pvalue,
                       intercept=NA,
                       intercept_se=NA,
                       inter_CIlower=NA,
                       inter_CIupper=NA,
                       intercept_pval=NA,
                       hetero_Q=result$Heter.Stat[1],
                       hetero_pvale=result$Heter.Stat[2])}
  else if (nrow(dat)>2 & all(dat$SNP %in% exp1$SNP)){rho <- ld_matrix(dat$SNP, pop = "EUR")
  result <- MendelianRandomization::mr_egger(mr_input(bx = dat$beta.exposure,
                                                      bxse = dat$se.exposure,
                                                      by = dat$beta.outcome,
                                                      byse = dat$se.outcome,
                                                      cor = rho))
  result <- data.frame(id=exp[i],
                       method="MR-Egger",
                       nsnp=result$SNPs,
                       b=result$Estimate,
                       se=result$StdError.Est,
                       CIlower=result$CILower.Est,
                       CIupper=result$CIUpper.Est,
                       pval=result$Pvalue.Est,
                       intercept=result$Intercept,
                       intercept_se=result$StdError.Int,
                       inter_CIlower=result$CILower.Int,
                       inter_CIupper=result$CIUpper.Int,
                       intercept_pval=result$Pvalue.Int,
                       hetero_Q=result$Heter.Stat[1],
                       hetero_pvale=result$Heter.Stat[2])
  result2 <- MendelianRandomization::mr_ivw(mr_input(bx = dat$beta.exposure,
                                                     bxse = dat$se.exposure,
                                                     by = dat$beta.outcome,
                                                     byse = dat$se.outcome,
                                                     cor = rho))
  result2 <- data.frame(id=exp[i],
                        method="IVW",
                        nsnp=result2$SNPs,
                        b=result2$Estimate,
                        se=result2$StdError,
                        CIlower=result2$CILower,
                        CIupper=result2$CIUpper,
                        pval=result2$Pvalue,
                        intercept=NA,
                        intercept_se=NA,
                        inter_CIlower=NA,
                        inter_CIupper=NA,
                        intercept_pval=NA,
                        hetero_Q=result2$Heter.Stat[1],
                        hetero_pvale=result2$Heter.Stat[2])
  result <- rbind(result,result2)}
  else if (nrow(dat)>2 & any(dat$SNP %in% exp2$SNP)){result <- MendelianRandomization::mr_egger(mr_input(bx = dat$beta.exposure,
                                                                                                         bxse = dat$se.exposure,
                                                                                                         by = dat$beta.outcome,
                                                                                                         byse = dat$se.outcome))
  result <- data.frame(id=exp[i],
                       method="MR-Egger",
                       nsnp=result$SNPs,
                       b=result$Estimate,
                       se=result$StdError.Est,
                       CIlower=result$CILower.Est,
                       CIupper=result$CIUpper.Est,
                       pval=result$Pvalue.Est,
                       intercept=result$Intercept,
                       intercept_se=result$StdError.Int,
                       inter_CIlower=result$CILower.Int,
                       inter_CIupper=result$CIUpper.Int,
                       intercept_pval=result$Pvalue.Int,
                       hetero_Q=result$Heter.Stat[1],
                       hetero_pvale=result$Heter.Stat[2])
  result2 <- MendelianRandomization::mr_ivw(mr_input(bx = dat$beta.exposure,
                                                     bxse = dat$se.exposure,
                                                     by = dat$beta.outcome,
                                                     byse = dat$se.outcome))
  result2 <- data.frame(id=exp[i],
                        method="IVW",
                        nsnp=result2$SNPs,
                        b=result2$Estimate,
                        se=result2$StdError,
                        CIlower=result2$CILower,
                        CIupper=result2$CIUpper,
                        pval=result2$Pvalue,
                        intercept=NA,
                        intercept_se=NA,
                        inter_CIlower=NA,
                        inter_CIupper=NA,
                        intercept_pval=NA,
                        hetero_Q=result2$Heter.Stat[1],
                        hetero_pvale=result2$Heter.Stat[2])
  result <- rbind(result,result2)}
  mrres <- rbind(mrres,result)
}
write.table(mrres, file="/Users/tq20202/Desktop/project2/asthma/ea/mrres_ea.txt", row.names=F, col.names=T, sep="\t", quote=F)


######African ancestry###########
rm(list=ls())
gc()

library(readxl)
library(TwoSampleMR)
library(MendelianRandomization)

exposure_samplesize <- 1871
outcome_samplesize <- 32653
ncase <- 5054
ncontrol <- 27599
prevalence <- ncase/outcome_samplesize

aa <- read_excel("/Users/tq20202/Desktop/project2/asthma/aa/conditional_aa.xlsx")

fname<-list.files("/Users/tq20202/Desktop/project2/aa_summary")
id=strsplit(fname[1],"[.]")[[1]][2]
name<-cbind(filename=fname[1],id)
for ( i in 2: length(fname)){
  temp<-cbind(fname[i],strsplit(fname[i],"[.]")[[1]][2])
  name<-rbind(name,temp)
}
name<-as.data.frame(name)

setwd("/Users/tq20202/Desktop/project2/aa_summary")
aares<-c()
for (i in 1: nrow(aa)){
  f<-name[name$id==as.character(aa[i,1]),1]
  if(length(f)==0){pi=rep(as.character(aa[i,1]),14)}
  else{temp <- read.table(file=f, header=T, sep="\t", quote="", stringsAsFactors=F, check.names=F)
  temp<-rbind(colnames(temp),temp[1:nrow(temp),])
  colnames(temp)<-c("CHROM","POS","ID",	"REF","ALT","A1",	"A1_FREQ","TEST",	"OBS_CT",	"BETA",	"SE",	"T_STAT",	"P",	"ERRCODE")
  pi<-temp[temp$ID==aa[i,10]$`Conditional independent cis-SNP`,]
  pi<-pi[1,]}
  aares<-rbind(aares,pi)
}
write.table(aares, file="/Users/tq20202/Desktop/project2/asthma/aa/aares.txt", row.names=F, col.names=T, sep="\t", quote=F)

aares<-cbind(aa[,1:5],aares)
aapqtl<-aares[aares$TEST=="ADD",] ##match summary stats
aapqtl<-cbind(aapqtl[,1:12],factor=0,ea=0,oa=0,aapqtl[,13:18])
aapqtl$A1_FREQ<-as.numeric(aapqtl$A1_FREQ)
for (i in 1:nrow(aapqtl)){
  ref<-aapqtl[i,"REF"]
  alt<-aapqtl[i,"ALT"]
  if (aapqtl[i,"A1"]==ref){aapqtl[i,"factor"]=-1;aapqtl[i,"ea"]=aapqtl[i,"A1"];aapqtl[i,"oa"]=aapqtl[i,"ALT"]}
  if (aapqtl[i,"A1"]==alt){aapqtl[i,"factor"]=1;aapqtl[i,"ea"]=aapqtl[i,"A1"];aapqtl[i,"oa"]=aapqtl[i,"REF"]}
}
unique(aapqtl$factor)

###format data
tmp <- data.frame(
  SNP = aapqtl$ID,
  beta = aapqtl$BETA,
  se = aapqtl$SE,
  effect_allele = aapqtl$ea,
  eaf = aapqtl$A1_FREQ,
  other_allele = aapqtl$oa,
  phenotype_id = aapqtl$`SOMAmer ID`,
  pval = aapqtl$P,
  chr = aapqtl$CHROM,
  pos = aapqtl$POS,
  samplesize = exposure_samplesize
)
exposure_dat <- format_data(tmp, type="exposure", effect_allele_col = "effect_allele",other_allele_col = "other_allele",phenotype_col = "phenotype_id",eaf_col = "eaf",samplesize_col = "samplesize",chr_col = "chr",pos_col = "pos")
#LD clumping
ld_reflookup <- function(rsid, pop='AFR')
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
write.table(exp1, file="/Users/tq20202/Desktop/project2/asthma/aa/exp1.txt", row.names=F, col.names=T, sep="\t", quote=F)
write.table(exp2, file="/Users/tq20202/Desktop/project2/asthma/aa/exp2.txt", row.names=F, col.names=T, sep="\t", quote=F)
exp1 <-clump_data(exp1,clump_kb=10000,clump_r2=0.6,clump_p1=1,clump_p2=1,pop="AFR")
exposure_dat<-rbind(exp1,exp2)

###F-statistics for each SNPs
exposure_dat <- cbind(exposure_dat,fstatistics=1)
for (s in 1:nrow(exposure_dat)){
  z <- exposure_dat[s,2]/exposure_dat[s,3]
  pve <- z^2/(exposure_samplesize+z^2)
  exposure_dat[s,"fstatistics"] <- (exposure_samplesize-2)*pve/(1-pve)
}
print(min(exposure_dat$fstatistics))
print(max(exposure_dat$fstatistics))
write.table(exposure_dat, file="/Users/tq20202/Desktop/project2/asthma/aa/exposure_dat_aa.txt", row.names=F, col.names=T, sep="\t", quote=F)

exposure_dat <- exposure_dat[exposure_dat$fstatistics>10,]
length(unique(exposure_dat$SNP))
length(unique(exposure_dat$exposure))

##outcome
Asthma_Bothsex_afr_inv_var_meta_GBMI_052021 <- read.delim("/Users/tq20202/outcomedata/Asthma_Bothsex_afr_inv_var_meta_GBMI_052021.txt", header=FALSE, comment.char="#")
outdate <- Asthma_Bothsex_afr_inv_var_meta_GBMI_052021[Asthma_Bothsex_afr_inv_var_meta_GBMI_052021$V5 %in% exposure_dat$SNP,]
length(unique(outdate$V5))
nosnp<-setdiff(exposure_dat$SNP,outdate$V5)
mexp<-exposure_dat[exposure_dat$SNP %in% nosnp,c("SNP","effect_allele.exposure","eaf.exposure","other_allele.exposure","beta.exposure")]
mexp<-mexp[!duplicated(mexp$SNP),]
######obtain unmatched snp outcome with proxy snp
aaproxy <- read.table(file="/Users/tq20202/Desktop/project2/AFR_pQTL_proxy.txt", header=T, sep="\t", quote="", stringsAsFactors=F, check.names=F)
aaproxy<-aaproxy[aaproxy$snp!=aaproxy$proxy,]
aaproxy<-aaproxy[aaproxy$snp %in% nosnp,]
s<-unique(aaproxy$snp)
proxydata<-c()
for (i in 1:length(s)){
  snp=s[i]
  d<-c()
  d2<-c()
  me<-mexp[mexp$SNP==snp,"effect_allele.exposure"]
  mef<-mexp[mexp$SNP==snp,"eaf.exposure"]
  mo<-mexp[mexp$SNP==snp,"other_allele.exposure"]
  temp<-aaproxy[aaproxy$snp==snp,"proxy"]
  out<-Asthma_Bothsex_afr_inv_var_meta_GBMI_052021[Asthma_Bothsex_afr_inv_var_meta_GBMI_052021$V5 %in% temp,]
  if (nrow(out)==0){
    d<-rep(NA,18)
  } else if (me %in% out$V4 & min(abs(out[out$V4==me,"V6"]-mef))<0.1){
    temp<-out[out$V4==me,]
    minf=min(abs(temp$V6-mef))
    d<-temp[abs(temp$V6-mef)==minf,]
    d<-cbind(d,snp)
  } else if (mo %in% out$V4 & min(abs((1-out[out$V4==mo,"V6"])-mef))<0.1){
    temp<-out[out$V4==mo,]
    minf=min(abs((1-temp$V6)-mef))
    d2<-temp[abs((1-temp$V6)-mef)==minf,]
    d2$V7=-d2$V7
    d2<-cbind(d2,snp)
  } else {
    minf=min(abs(out$V6-mef))
    d<-out[abs(out$V6-mef)==minf,]
    d<-cbind(d,snp)
  }
  proxydata<-rbind(proxydata,d,d2)
}
proxy<-na.omit(proxydata[,c(1:9,18)])
msnp<-c()
for (i in 1:nrow(proxy)){
  temp=proxy[i,10]
  temp2<-mexp[mexp$SNP==temp,2:5]
  msnp<-rbind(msnp,temp2)
}
proxy<-cbind(proxy,msnp)

proxy<-proxy[abs(proxy$V6-proxy$eaf.exposure)<0.1,]
length(unique(proxy$snp))
write.table(proxy, file="/Users/tq20202/Desktop/project2/asthma/aa/proxy_aa.txt", row.names=F, col.names=T, sep="\t", quote=F)

proxy<-proxy[!duplicated(proxy$snp),]
proxy_data<-data.frame(SNP=proxy$snp,
                       effect_allele=proxy$effect_allele.exposure,
                       eaf=proxy$eaf.exposure,
                       other_allele=proxy$other_allele.exposure,
                       beta=proxy$V7,
                       se=proxy$V8,
                       pval=proxy$V9)
out_data<-data.frame(SNP=outdate$V5,
                     effect_allele=outdate$V4,
                     eaf=outdate$V6,
                     other_allele=outdate$V3,
                     beta=outdate$V7,
                     se=outdate$V8,
                     pval=outdate$V9)
outcome_dat <- rbind(out_data,proxy_data)
length(unique(outcome_dat$SNP))
outcome_dat  <- format_data(outcome_dat , type="outcome")
write.table(outcome_dat, file="/Users/tq20202/Desktop/project2/asthma/aa/outcome_dat_aa.txt", row.names=F, col.names=T, sep="\t", quote=F)

#####2SMR for matched data#################################################
inter<-intersect(exposure_dat$SNP,outcome_dat$SNP)
exposure_dat <- exposure_dat[exposure_dat$SNP %in% inter,]
outcome_dat <- outcome_dat[outcome_dat$SNP %in% inter,]
harmdat <- harmonise_data(exposure_dat, outcome_dat, action=2)
sum(abs(harmdat$eaf.exposure-harmdat$eaf.outcome)<0.1)
tharm=harmdat
tharm$test = tharm$effect_allele.exposure == tharm$effect_allele.outcome
nrow(tharm[tharm$test==FALSE,])
harmdat$mr_keep<-TRUE  #consider ambiguous palindromic SNPs (for the "MendelianRandomization" package)

test<-harmdat
test$new<- paste(test$SNP,test$exposure,sep="_")
test<-test[,-24]
colnames(test)[31]<-"exposure"
test<-test[,-26]
test <- format_data(test, type="exposure", snp_col = "SNP", pval_col = "pval.exposure",
                    beta_col = "beta.exposure", se_col = "se.exposure",
                    effect_allele_col = "effect_allele.exposure",
                    other_allele_col = "other_allele.exposure",
                    phenotype_col = "exposure",eaf_col = "eaf.exposure",    
                    samplesize_col = "samplesize.exposure")
testout <- rbind(out_data,proxy_data)
testout  <- format_data(testout , type="outcome")
inter<-intersect(test$SNP,testout$SNP)
test <- test[test$SNP %in% inter,]
testout <- testout[testout$SNP %in% inter,]
testharm <- harmonise_data(test, testout, action=2)
testharm$mr_keep<-TRUE
testres<-mr(testharm)
###steiger filtering
testres <- generate_odds_ratios(testres)
testres <- cbind(SNP=1,testres,beta.exp=1,beta.out=1,se.exp=1,se.out=1,eaf.exp=1,eaf.out=1,p.exp=1,p.out=1)
for (r in 1: nrow(testharm)){
  testres[r,1] <- testharm[testharm$id.exposure==testres[r,"id.exposure"],"SNP"]
  testres[r,16:23] <- testharm[testharm$id.exposure==testres[r,"id.exposure"],c("beta.exposure","beta.outcome","se.exposure",
                                                                                "se.outcome","eaf.exposure","eaf.outcome",
                                                                                "pval.exposure","pval.outcome")]
}
testres <- na.omit(testres)
testres <- cbind(testres,rsq.exposure=1,rsq.outcome=1,samplesize.exposure=exposure_samplesize,samplesize.outcome=outcome_samplesize)
testres$rsq.exposure <- (get_r_from_pn(p=testres$p.exp, n=exposure_samplesize))^2
testres$rsq.outcome <- (get_r_from_lor(lor=log(testres$or),af=testres$eaf.out,ncase=ncase,ncontrol=ncontrol,prevalence=prevalence))^2
st <- psych::r.test(
  n = testres$samplesize.exposure, 
  n2 = testres$samplesize.outcome, 
  r12 = sqrt(testres$rsq.exposure), 
  r34 = sqrt(testres$rsq.outcome)
)
testres$steiger_dir <- testres$rsq.exposure > testres$rsq.outcome
testres$steiger_pval <- pnorm(-abs(st$z)) * 2

harmdat$match<-paste(harmdat$SNP,harmdat$exposure,sep="_")
harmdat$steiger_dir <- 0
harmdat$steiger_pval <- 0
for (i in 1:nrow(harmdat)){
  harmdat[i,33] <- testres[testres$exposure==harmdat[i,"match"],"steiger_dir"]
  harmdat[i,34] <- testres[testres$exposure==harmdat[i,"match"],"steiger_pval"]
}
harmdat<-harmdat[,-32]
harmdat<-harmdat[harmdat$steiger_dir==1,]
length(unique(harmdat$exposure))
write.table(harmdat, file="/Users/tq20202/Desktop/project2/asthma/aa/harmdat_aa.txt", row.names=F, col.names=T, sep="\t", quote=F)

#Generalized IVW  and MR-Egger for the "MendelianRandomization" package
exp<-unique(harmdat$exposure)
mrres <- c()
for (i in 1:length(exp)){
  dat <- harmdat[harmdat$exposure==exp[i],]
  if (nrow(dat)==1){
    result <- mr_wald_ratio(b_exp= dat$beta.exposure, b_out=dat$beta.outcome, se_exp= dat$se.exposure, se_out= dat$se.outcome)
    result <- data.frame(id=exp[i],method="Wald_ratio",nsnp=result$nsnp,b=result$b,se=result$se,CIlower=NA,CIupper=NA,
                         pval=result$pval,intercept=NA,intercept_se=NA,inter_CIlower=NA,inter_CIupper=NA,
                         intercept_pval=NA,hetero_Q=NA,hetero_pvale=NA)}
  else if (nrow(dat)==2 & all(dat$SNP %in% exp1$SNP)){rho <- ld_matrix(dat$SNP, pop = "AFR")
  result <- MendelianRandomization::mr_ivw(mr_input(bx = dat$beta.exposure,
                                                    bxse = dat$se.exposure,
                                                    by = dat$beta.outcome,
                                                    byse = dat$se.outcome,
                                                    cor = rho))
  result <- data.frame(id=exp[i],
                       method="IVW",
                       nsnp=result$SNPs,
                       b=result$Estimate,
                       se=result$StdError,
                       CIlower=result$CILower,
                       CIupper=result$CIUpper,
                       pval=result$Pvalue,
                       intercept=NA,
                       intercept_se=NA,
                       inter_CIlower=NA,
                       inter_CIupper=NA,
                       intercept_pval=NA,
                       hetero_Q=result$Heter.Stat[1],
                       hetero_pvale=result$Heter.Stat[2])}
  else if (nrow(dat)==2 & any(dat$SNP %in% exp2$SNP)){result <- MendelianRandomization::mr_ivw(mr_input(bx = dat$beta.exposure,
                                                                                                        bxse = dat$se.exposure,
                                                                                                        by = dat$beta.outcome,
                                                                                                        byse = dat$se.outcome))
  result <- data.frame(id=exp[i],
                       method="IVW",
                       nsnp=result$SNPs,
                       b=result$Estimate,
                       se=result$StdError,
                       CIlower=result$CILower,
                       CIupper=result$CIUpper,
                       pval=result$Pvalue,
                       intercept=NA,
                       intercept_se=NA,
                       inter_CIlower=NA,
                       inter_CIupper=NA,
                       intercept_pval=NA,
                       hetero_Q=result$Heter.Stat[1],
                       hetero_pvale=result$Heter.Stat[2])}
  else if (nrow(dat)>2 & all(dat$SNP %in% exp1$SNP)){rho <- ld_matrix(dat$SNP, pop = "AFR")
  result <- MendelianRandomization::mr_egger(mr_input(bx = dat$beta.exposure,
                                                      bxse = dat$se.exposure,
                                                      by = dat$beta.outcome,
                                                      byse = dat$se.outcome,
                                                      cor = rho))
  result <- data.frame(id=exp[i],
                       method="MR-Egger",
                       nsnp=result$SNPs,
                       b=result$Estimate,
                       se=result$StdError.Est,
                       CIlower=result$CILower.Est,
                       CIupper=result$CIUpper.Est,
                       pval=result$Pvalue.Est,
                       intercept=result$Intercept,
                       intercept_se=result$StdError.Int,
                       inter_CIlower=result$CILower.Int,
                       inter_CIupper=result$CIUpper.Int,
                       intercept_pval=result$Pvalue.Int,
                       hetero_Q=result$Heter.Stat[1],
                       hetero_pvale=result$Heter.Stat[2])
  result2 <- MendelianRandomization::mr_ivw(mr_input(bx = dat$beta.exposure,
                                                     bxse = dat$se.exposure,
                                                     by = dat$beta.outcome,
                                                     byse = dat$se.outcome,
                                                     cor = rho))
  result2 <- data.frame(id=exp[i],
                        method="IVW",
                        nsnp=result2$SNPs,
                        b=result2$Estimate,
                        se=result2$StdError,
                        CIlower=result2$CILower,
                        CIupper=result2$CIUpper,
                        pval=result2$Pvalue,
                        intercept=NA,
                        intercept_se=NA,
                        inter_CIlower=NA,
                        inter_CIupper=NA,
                        intercept_pval=NA,
                        hetero_Q=result2$Heter.Stat[1],
                        hetero_pvale=result2$Heter.Stat[2])
  result <- rbind(result,result2)}
  else if (nrow(dat)>2 & any(dat$SNP %in% exp2$SNP)){result <- MendelianRandomization::mr_egger(mr_input(bx = dat$beta.exposure,
                                                                                                         bxse = dat$se.exposure,
                                                                                                         by = dat$beta.outcome,
                                                                                                         byse = dat$se.outcome))
  result <- data.frame(id=exp[i],
                       method="MR-Egger",
                       nsnp=result$SNPs,
                       b=result$Estimate,
                       se=result$StdError.Est,
                       CIlower=result$CILower.Est,
                       CIupper=result$CIUpper.Est,
                       pval=result$Pvalue.Est,
                       intercept=result$Intercept,
                       intercept_se=result$StdError.Int,
                       inter_CIlower=result$CILower.Int,
                       inter_CIupper=result$CIUpper.Int,
                       intercept_pval=result$Pvalue.Int,
                       hetero_Q=result$Heter.Stat[1],
                       hetero_pvale=result$Heter.Stat[2])
  result2 <- MendelianRandomization::mr_ivw(mr_input(bx = dat$beta.exposure,
                                                     bxse = dat$se.exposure,
                                                     by = dat$beta.outcome,
                                                     byse = dat$se.outcome))
  result2 <- data.frame(id=exp[i],
                        method="IVW",
                        nsnp=result2$SNPs,
                        b=result2$Estimate,
                        se=result2$StdError,
                        CIlower=result2$CILower,
                        CIupper=result2$CIUpper,
                        pval=result2$Pvalue,
                        intercept=NA,
                        intercept_se=NA,
                        inter_CIlower=NA,
                        inter_CIupper=NA,
                        intercept_pval=NA,
                        hetero_Q=result2$Heter.Stat[1],
                        hetero_pvale=result2$Heter.Stat[2])
  result <- rbind(result,result2)}
  mrres <- rbind(mrres,result)
}
write.table(mrres, file="/Users/tq20202/Desktop/project2/asthma/aa/mrres_aa.txt", row.names=F, col.names=T, sep="\t", quote=F)

#######LDcheck and colocalization for sensitivity analysis###################
####LDcheck#######
##outcome
Asthma_Bothsex_nfe_inv_var_meta_GBMI_052021 <- read.delim("/Users/tq20202/outcomedata/Asthma_Bothsex_nfe_inv_var_meta_GBMI_052021.txt", header=FALSE, comment.char="#")
##exposure
harmdat_dat <- read.table(file="/Users/tq20202/Desktop/project2/asthma/ea/harmdat_ea.txt", header=T, sep="\t", quote="", stringsAsFactors=F, check.names=F)
harmdat_dat$chr <- harmdat_dat$chr.exposure
harmdat_dat$pos1 <- harmdat_dat$pos.exposure-500000
harmdat_dat$pos2 <- harmdat_dat$pos.exposure+500000
res<-read.table(file="/Users/tq20202/Desktop/project2/asthma/ldcheck/mrres_ea.txt", header=T, sep="\t", quote="", stringsAsFactors=F, check.names=F)
res<-res[res$significant==T,]
harm<-harmdat_dat[harmdat_dat$exposure %in% res$id,]
snpfile <- harm
#####find genetic association info of all snps within window########
#####get LD_Matrix and LD_r2########################################
ldres <- c()
for (i in 1:nrow(snpfile)){
  print(i)
  rsid <- snpfile$SNP[i]
  c<-snpfile$chr[i]
  p1<-snpfile$pos1[i]
  p2<-snpfile$pos2[i]
  assoc<-Asthma_Bothsex_nfe_inv_var_meta_GBMI_052021[Asthma_Bothsex_nfe_inv_var_meta_GBMI_052021$V1==c &Asthma_Bothsex_nfe_inv_var_meta_GBMI_052021$V2>p1 & Asthma_Bothsex_nfe_inv_var_meta_GBMI_052021$V2<p2,]
  if (nrow(assoc)!=0){
    assoc <- assoc[order(assoc$V9),]
    data <- assoc[assoc$V9<1E-3,]
    if(nrow(data)>=500){data <- data[1:499,] }
    print(nrow(data))
    
    if (nrow(data)!=0){
      data<-na.omit(data[,1:9])
      snp <- append(data$V5, rsid)
      a <- NULL
      attempts <- 0
      while(attempts<=10){    
        a <- ld_matrix(variants=snp,pop = "EUR")
        if(is.null(a)){attempts<-attempts+1}else{break}
      }
      if(is.null(nrow(a))==TRUE){c<-cbind(snp=rsid, ld_snp=rsid, ld_r2=1)} 
      else {col.index <- which(grepl(rsid,colnames(a)))
      if (length(col.index)>0){
        b <- (a[,col.index])^2
        b <- b[order(b)]
        b <- b[(length(b)-1)] 
        c <- cbind(snp=rsid, ld_snp=names(b), ld_r2=as.numeric(b))
      } else {c <- cbind(snp=rsid, ld_snp="NA", ld_r2="NA")}
      }
    } else {
      c <- c(snp=rsid, ld_snp="NA",ld_r2="NA")}
  } else {c <- c(snp=rsid, ld_snp="NA",ld_r2="NA")}
  
  ldres <- rbind(ldres,c)
}
sum(harm$SNP==ldres[,1])
r<-cbind(harm,ldres[,2:3])
write.table(r, file="/Users/tq20202/Desktop/project2/asthma/ldcheck/LDcheck_asthma_ea.txt", row.names=F, col.names=T, sep="\t", quote=F)

e<-unique(ld$exposure)
result<-c()
for (i in 1: length(e)){
  p<-e[i]
  temp<-ld[ld$exposure==p,]
  if (sum(temp$ldcheck)==0){l=FALSE}else{l=TRUE}
  res<-cbind(p,l)
  result<-rbind(result,res)
}
write.table(result, file="/Users/tq20202/Desktop/project2/asthma/ldcheck/LDres_asthma_ea.txt", row.names=F, col.names=T, sep="\t", quote=F)

##outcome
Asthma_Bothsex_afr_inv_var_meta_GBMI_052021 <- read.delim("/Users/tq20202/outcomedata/Asthma_Bothsex_afr_inv_var_meta_GBMI_052021.txt", header=FALSE, comment.char="#")
##exposure
harmdat_dat <- read.table(file="/Users/tq20202/Desktop/project2/asthma/aa/harmdat_aa.txt", header=T, sep="\t", quote="", stringsAsFactors=F, check.names=F)
harmdat_dat$chr <- harmdat_dat$chr.exposure
harmdat_dat$pos1 <- harmdat_dat$pos.exposure-500000
harmdat_dat$pos2 <- harmdat_dat$pos.exposure+500000
res<-read.table(file="/Users/tq20202/Desktop/project2/asthma/ldcheck/mrres_aa.txt", header=T, sep="\t", quote="", stringsAsFactors=F, check.names=F)
res<-res[res$significant==T,]
harm<-harmdat_dat[harmdat_dat$exposure %in% res$id,]
snpfile <- harm
#####find genetic association info of all snps within window########
#####get LD_Matrix and LD_r2########################################
ldres <- c()
for (i in 1:nrow(snpfile)){
  print(i)
  rsid <- snpfile$SNP[i]
  c<-snpfile$chr[i]
  p1<-snpfile$pos1[i]
  p2<-snpfile$pos2[i]
  assoc<-Asthma_Bothsex_afr_inv_var_meta_GBMI_052021[Asthma_Bothsex_afr_inv_var_meta_GBMI_052021$V1==c &Asthma_Bothsex_afr_inv_var_meta_GBMI_052021$V2>p1 & Asthma_Bothsex_afr_inv_var_meta_GBMI_052021$V2<p2,]
  if (nrow(assoc)!=0){
    assoc <- assoc[order(assoc$V9),]
    data <- assoc[assoc$V9<1E-3,]
    if(nrow(data)>=500){data <- data[1:499,] }
    print(nrow(data))
    
    if (nrow(data)!=0){
      data<-na.omit(data[,1:9])
      snp <- append(data$V5, rsid)
      a <- NULL
      attempts <- 0
      while(attempts<=10){    
        a <- ld_matrix(variants=snp,pop = "AFR")
        if(is.null(a)){attempts<-attempts+1}else{break}
      }
      if(is.null(nrow(a))==TRUE){c<-cbind(snp=rsid, ld_snp=rsid, ld_r2=1)} 
      else {col.index <- which(grepl(rsid,colnames(a)))
      if (length(col.index)>0){
        b <- (a[,col.index])^2
        b <- b[order(b)]
        b <- b[(length(b)-1)] 
        c <- cbind(snp=rsid, ld_snp=names(b), ld_r2=as.numeric(b))
      } else {c <- cbind(snp=rsid, ld_snp="NA", ld_r2="NA")}
      }
    } else {
      c <- c(snp=rsid, ld_snp="NA",ld_r2="NA")}
  } else {c <- c(snp=rsid, ld_snp="NA",ld_r2="NA")}
  
  ldres <- rbind(ldres,c)
}
sum(harm$SNP==ldres[,1])
r<-cbind(harm,ldres[,2:3])
write.table(r, file="/Users/tq20202/Desktop/project2/asthma/ldcheck/LDcheck_asthma_aa.txt", row.names=F, col.names=T, sep="\t", quote=F)

e<-unique(ld$exposure)
result<-c()
for (i in 1: length(e)){
  p<-e[i]
  temp<-ld[ld$exposure==p,]
  if (sum(temp$ldcheck)==0){l=FALSE}else{l=TRUE}
  res<-cbind(p,l)
  result<-rbind(result,res)
}
write.table(result, file="/Users/tq20202/Desktop/project2/asthma/ldcheck/LDres_asthma_ea.txt", row.names=F, col.names=T, sep="\t", quote=F)

####colocalization#####
###use pwcoco######
sigea <- read_excel("sigres_ea.xlsx")
sigaa <- read_excel("sigres_aa.xlsx")

s1<-"pwcoco --bfile /newhome/tq20202/pwcoco/pwcoco-master/EUR/chr"
s2<-"/chr"
s3<-'  --sum_stats1 "prepdata/'
s4<-'_ea.txt" --sum_stats2 "prepdata/'
s5<-'_ea.txt" --out "'
s6<-'" --out_cond'

eascript<-c()
for (i in 1:nrow(sigea)){
  chr<-sigea[i,6]$chromosome
  g<-sigea[i,5]$`Entrez Gene Symbol`
  d<-sigea[i,7]$disease
  temp<-paste(s1,chr,s2,chr,s3,g,s4,g,"_",d,s5,"outcome/",g,"_ea_",d,"_res",s6,sep="")
  eascript <- rbind(eascript,temp)
}


s1<-"pwcoco --bfile /newhome/tq20202/pwcoco/pwcoco-master/AFR/chr"
s2<-"/chr"
s3<-'  --sum_stats1 "prepdata/'
s4<-'_aa.txt" --sum_stats2 "prepdata/'
s5<-'_aa.txt" --out "'
s6<-'" --out_cond'

aascript<-c()
for (i in 1:nrow(sigaa)){
  chr<-sigaa[i,6]$chromosome
  g<-sigaa[i,5]$`Entrez Gene Symbol`
  d<-sigaa[i,7]$disease
  efile<-paste(g,"_aa.txt",sep="")
  ofile<-paste(g,"_",d,"_aa.txt",sep="")
  temp<-paste(s1,chr,s2,chr,s3,g,s4,g,"_",d,s5,"outcome/",g,"_aa_",d,"_res",s6,sep="")
  aascript <- rbind(aascript,temp)
}

write.table(eascript, file="eascript.txt", row.names=F, col.names=T, sep="\t", quote=F)
write.table(aascript, file="aascript.txt", row.names=F, col.names=T, sep="\t", quote=F)

f<-list.files("out")
setwd("/Users/tq20202/Desktop/project2/pwcoco/out")
pres<-c()
for (i in 1:length(f)){
  g<-strsplit(f[i],"_")[[1]][1]
  an<-strsplit(f[i],"_")[[1]][2]
  d<-strsplit(f[i],"_")[[1]][3]
  temp<-data.frame(gene=g,ancestry=an,disease=d,H4=0,pwcoco=FALSE)
  data<-read.table(file=f[i], header=T, sep="\t", quote="", stringsAsFactors=F, check.names=F)
  temp$H4<-max(data$H4)
  if(temp$H4>0.7){temp$pwcoco=TRUE}
  pres<-rbind(pres,temp)
}
write.table(pres, file="/Users/tq20202/Desktop/project2/pwcoco/pwcoco_res.txt", row.names=F, col.names=T, sep="\t", quote=F)









