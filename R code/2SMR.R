library(TwoSampleMR)
library(MendelianRandomization)
library(ieugwasr)

#####2SMR for matched data#################################################
inter<-intersect(exposure_dat$SNP,outcome_dat$SNP)
exposure_dat <- exposure_dat[exposure_dat$SNP %in% inter,]
outcome_dat <- outcome_dat[outcome_dat$SNP %in% inter,]
harmdat <- harmonise_data(exposure_dat, outcome_dat, action=2)
harmdat$mr_keep<-TRUE  #consider ambiguous palindromic SNPs (for the "MendelianRandomization" package)

###steiger filtering
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

#####heterogenity and pleiotropy test summary##############
mr1<-mrres[mrres$method=="Wald_ratio",]
mr1$heterogenity<-NA
mr1$pleiotropy<-NA
mr2<-mrres[mrres$method=="IVW" & mrres$nsnp==2,]
for (i in 1:nrow(mr2)){
  if (mr2[i,"hetero_pvale"]<0.05){mr2[i,"heterogenity"]<-"TRUE"}else{mr2[i,"heterogenity"]<-"FALSE"}
}
mr2$pleiotropy<-NA
mr3<-mrres[mrres$method=="IVW" & mrres$nsnp>2,]
for (i in 1:nrow(mr3)){
  if (mr3[i,"hetero_pvale"]<0.05){mr3[i,"heterogenity"]<-"TRUE"}else{mr3[i,"heterogenity"]<-"FALSE"}
}
mr4<-mrres[mrres$method=="MR-Egger",]
sum(mr3$id==mr4$id)==nrow(mr3)
for (i in 1:nrow(mr4)){
  if (mr4[i,"intercept_pval"]>0.05){mr3[i,"pleiotropy"]<-"FALSE"}else{mr3[i,"pleiotropy"]<-"TRUE"}
}
mr<-rbind(mr1,mr2,mr3)

#####FDR Adjustment for Pvalue######
mr$fdr <- p.adjust(mr$pval, method = "fdr", n = length(mr$pval))




