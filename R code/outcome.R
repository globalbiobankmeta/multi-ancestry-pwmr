library(TwoSampleMR)

##outcome
Asthma_Bothsex_nfe_inv_var_meta_GBMI_052021 <- read.delim("Asthma_Bothsex_nfe_inv_var_meta_GBMI_052021.txt", header=FALSE, comment.char="#")
outdate <- Asthma_Bothsex_nfe_inv_var_meta_GBMI_052021[Asthma_Bothsex_nfe_inv_var_meta_GBMI_052021$V5 %in% exposure_dat$SNP,]

######obtain unmatched snp outcome with proxy snp
nosnp<-setdiff(exposure_dat$SNP,outdate$V5)
mexp<-exposure_dat[exposure_dat$SNP %in% nosnp,c("SNP","effect_allele.exposure","eaf.exposure","other_allele.exposure","beta.exposure")]
mexp<-mexp[!duplicated(mexp$SNP),]
eaproxy <- read.table(file="EUR_pQTL_proxy.txt", header=T, sep="\t", quote="", stringsAsFactors=F, check.names=F)
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
outcome_dat  <- format_data(outcome_dat , type="outcome")
