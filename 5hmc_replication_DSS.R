library(DSS)
library(bsseq)

######description####
#The group of samples in here was processed with Nugen kit 

#Oxidatiaon was done before bisulfite conversion(called hm samples here)
#Yield methly calls from Bismark from these samples contains both 5hmc and 5mc 

#Comparison with replicated in either wildtype or knockout group will give us a panel of reproducible
#calls in either group


#####data read in and build BS object####
# Uhrf1 Knockout rep1(file very large)
data_MT1_hm <- read.table("data/MT1_hm_mc_hmc_bisulfite_trimmed_bismark_bt2.CpG_report_for_dss.txt")
colnames(data_MT1_hm)<- c("chr", "pos", "str", "N", "X")
# Uhrf1 Knockout rep2
data_MT2_hm <- read.table("data/MT2_hm_mc_hmc_bisulfite_trimmed_bismark_bt2.CpG_report_for_dss.txt")
colnames(data_MT2_hm)<- c("chr", "pos", "str", "N", "X")
# Uhrf1 wildtype rep1
data_WT1_hm <- read.table("data/WT1_hm_mc_hmc_bisulfite_trimmed_bismark_bt2.CpG_report_for_dss.txt")
colnames(data_WT1_hm)<- c("chr", "pos", "str", "N", "X")
# Uhrf1 wildtype rep2
data_WT2_hm <- read.table("data/WT2_hm_mc_hmc_bisulfite_trimmed_bismark_bt2.CpG_report_for_dss.txt")
colnames(data_WT2_hm)<- c("chr", "pos", "str", "N", "X")
# make BSseq objects 
BSobjMT <- makeBSseqData(list(data_MT1_hm, data_MT2_hm), c("MT1","MT2"))
BSobjWT <- makeBSseqData(list(data_WT1_hm, data_WT2_hm), c("WT1","WT2"))



###For whole-genome BS-seq data, perform DML test with smoothing ####
dmlTest_MT <- DMLtest(BSobjMT, group1=c("MT1"), group2=c("MT2"), smoothing=TRUE, smoothing.span=500)
dmlTest_WT <- DMLtest(BSobjWT, group1=c("WT1"), group2=c("WT2"), smoothing=TRUE, smoothing.span=500)
## call DML
# dmls <- callDML(dmlTest)
## call DML with a threshold
# dmls2 <- callDML(dmlTest, delta=0.2)




#### filter and save the data####
dmlTest_MT_filter=subset(dmlTest_MT,select=-c(diff.se,phi1,phi2))
dmlTest_WT_filter=subset(dmlTest_WT,select=-c(diff.se,phi1,phi2))

##adding average of mu0 and mu1
dmlTest_MT_filter$average<-(dmlTest_MT_filter$mu1+dmlTest_MT_filter$mu2)/2
dmlTest_WT_filter$average<-(dmlTest_WT_filter$mu1+dmlTest_WT_filter$mu2)/2

##filter out average value <0.05, filtering out noise 
dmlTest_MT_filter=subset(dmlTest_MT_filter,dmlTest_MT_filter$average>0.05)
dmlTest_MT_filter=subset(dmlTest_MT_filter,dmlTest_MT_filter$average>0.05)

##save the readout 
write.csv(dmlTest_MT_filter,file = "dmlTest_MT_filter.csv")
write.csv(dmlTest_WT_filter,file = "dmlTest_WT_filter.csv")


