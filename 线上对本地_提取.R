library(data.table)
library(TwoSampleMR)
GERD<-fread(file = "finngen_R9_M13_SPONDYLOLISTHESIS.gz",header = T)
head(GERD)
setwd("D:/桌面/腰椎滑脱/")
getwd()
#ebi-a-GCST90012877#阿兹海默
EEEEEEE<-extract_instruments(
  outcomes="ukb-b-10011",
  p1 = 5e-08,
  clump=TRUE, 
  r2=0.001,                  
  kb=10000,
  access_token= NULL) 
write.csv(EEEEEEE,file="EEEEEEE.csv")

GERD1<-merge(EEEEEEE,
             GERD,
            by.x = "SNP",
            by.y = "rsids")


write.csv(GERD1,file="GERD1.csv")

GERD1<-read_outcome_data(filename = "GERD1.csv",
                        sep = ",",
                        snp_col = "SNP",
                        beta_col = "beta",
                        se_col = "sebeta",
                        effect_allele_col = "alt",
                        other_allele_col = "ref",eaf_col = "af_alt",
                        pval_col = "pval")




mydata<-harmonise_data(
  exposure_dat=EEEEEEE,
  outcome_dat=GERD1,
  action= 2  
)                                #同方向纠正 
mydata<-subset(mydata,mr_keep)

################################################################################
#修改暴露和结局的列名称，为作图准备
mydata$outcome <- "Spondylolisthesis"
mydata$samplesize.outcome <-276730   #Spondylolisthesis/Spondylolysis

mydata$exposure <- "Years of schooling"
mydata$exposure <- "Job involves heavy manual or physical work"
mydata$exposure <- "Average total household income before tax"
mydata$exposure <- "Townsend deprivation index at recruitment"






################################################################################
res <-mr(mydata)
res  #结果计算

write.csv(mydata,file="mydata.csv")

