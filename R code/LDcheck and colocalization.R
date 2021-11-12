#######LDcheck and colocalization for sensitivity analysis###################
library(TwoSampleMR)
library(MendelianRandomization)
library(ieugwasr)

####LDcheck#######
##outcome
Asthma_Bothsex_nfe_inv_var_meta_GBMI_052021 <- read.delim("Asthma_Bothsex_nfe_inv_var_meta_GBMI_052021.txt", header=FALSE, comment.char="#")
##exposure
harmdat_dat <- read.table(file="harmdat_ea.txt", header=T, sep="\t", quote="", stringsAsFactors=F, check.names=F)
harmdat_dat$chr <- harmdat_dat$chr.exposure
harmdat_dat$pos1 <- harmdat_dat$pos.exposure-500000
harmdat_dat$pos2 <- harmdat_dat$pos.exposure+500000
res<-read.table(file="mrres_ea.txt", header=T, sep="\t", quote="", stringsAsFactors=F, check.names=F)
res<-res[res$fdr<0.05,]
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
ld<-cbind(harm,ldres[,2:3])

e<-unique(ld$exposure)
ldresult<-c()
for (i in 1: length(e)){
  p<-e[i]
  temp<-ld[ld$exposure==p,]
  if (sum(temp$ldcheck)==0){l=FALSE}else{l=TRUE}
  res<-cbind(p,l)
  ldresult<-rbind(ldresult,res)
}


####colocalization#####
###use pwcoco in linux######
if(pwcoco_out$H4>0.7){pwcoco_out$pwcoco=TRUE}else{pwcoco_out$pwcoco=FALSE}

