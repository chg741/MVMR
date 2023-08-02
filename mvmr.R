#读入数据包

library(TwoSampleMR)
library(data.table)

setwd("D:/桌面/腰椎滑脱/")


id_exposure <- c("ieu-a-1239","ukb-b-2002","ukb-b-7408","ukb-b-10011") 
#id_outcome <- "---"
exposure_dat <- mv_extract_exposures(id_exposure,clump_r2 = 0.001,clump_kb = 10000,pval_threshold = 5e-08,find_proxies = T)#三个显著相关的SNP提出
#dim(exposure_dat)
# #################################################################################
write.csv(exposure_dat, file="exposure_dat.csv")
exposure_dat<-read.csv("exposure_dat.csv", header = TRUE)

FINN_SL<-fread(file = "finngen_R9_M13_SPONDYLOLISTHESIS.gz",header = T)
outcome_dat<-merge(exposure_dat,FINN_SL,
      by.x = "SNP",by.y = "rsids")


outcome_dat<-clump_data(dat = outcome_dat,clump_kb = 10000,clump_r2 = 0.001,pop = "EUR")
head(outcome_dat)
outcome_dat<-outcome_dat[,-2]
outcome_dat<-outcome_dat[,c(1,12:13,15,17:19)]

#改列名
colnames(outcome_dat)[7] <- "eaf.outcome"
colnames(outcome_dat)[6] <- "se.outcome"
colnames(outcome_dat)[5] <- "beta.outcome"
colnames(outcome_dat)[4] <- "pval.outcome"
colnames(outcome_dat)[3] <- "effect_allele.outcome"
colnames(outcome_dat)[2] <- "other_allele.outcome"
head(outcome_dat)
#增加新列
outcome_dat$outcome<-"FINN"

outcome_dat$id.outcome<-"SL"

###################################################################################
#outcome_dat <- extract_outcome_data(exposure_dat$SNP, id_outcome,maf_threshold = 0.01,proxies = T)

mvdat <- mv_harmonise_data(exposure_dat, outcome_dat)
#################################################################################条件F值Conditional F-statistic
library(MVMR)

mv_input <- MVMR::format_mvmr(
  BXGs = mvdat$exposure_beta,
  seBXGs = mvdat$exposure_se,
  BYG = mvdat$outcome_beta,
  seBYG = mvdat$outcome_se,
  RSID = rownames(mvdat$exposure_beta)
)

F_strength <- strength_mvmr(mv_input)
################################################################################
res <- mv_multiple(mvdat,instrument_specific=T,pval_threshold=5E-8)
res
############################################
res_OR<-generate_odds_ratios(res$result)
res_OR
write.csv(res_OR, file="res_mvmr.csv")
##############################################
# res$result$lci95 <-res$result$b-1.96*res$result$se
# res$result$upi95 <-res$result$b+1.96*res$result$se
# res$result$or <-exp(res$result$b)
# res$result$or_lci95 <-exp(res$result$lci95)
# res$result$or_upi95 <-exp(res$result$upi95)
# res
# write.csv(res, file="res.csv")
#################################################################################
######LASOO挑选变量############################################################
mv_lasso_feature_selection(mvdat)
########################单个snp多变量的方法
mv_residual(
  mvdat,
  intercept = FALSE,
  instrument_specific = FALSE,
  pval_threshold = 5e-08,
  plots = FALSE
)

mv_multiple(
  mvdat,
  intercept = FALSE,
  instrument_specific = FALSE,
  pval_threshold = 5e-08,
  plots = FALSE
)

mv_subset(
  mvdat,
  features = mv_lasso_feature_selection(mvdat),
  intercept = FALSE,
  instrument_specific = FALSE,
  pval_threshold = 5e-08,
  plots = FALSE
)
