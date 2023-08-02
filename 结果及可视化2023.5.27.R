library(TwoSampleMR)
library(MRPRESSO)
setwd("D:\\桌面\\data all/INCOME/")
getwd()
# mr_method_list()
Mydata<-read.csv("mydata.csv",header = T)
Mydata<-subset(Mydata,mr_keep)
res <-mr(Mydata,method_list=c("mr_egger_regression","mr_weighted_median", "mr_ivw","mr_raps","mr_simple_mode","mr_weighted_mode"))

res  #结果计算

# ##########################################################################
#  异质性，敏感性，多效性分析
heterogeneity<- mr_heterogeneity(Mydata, method_list=c("mr_egger_regression", "mr_ivw"))  #异质性检验——IVWorMRegger   #I2=[Q-(K-1)]/Q
heterogeneity   


pleio <- mr_pleiotropy_test(Mydata)
pleio   #多效性检验——MR egger

####################################################################################
k=nrow(Mydata)

if (((k/0.05)+100) > 1000) {
  Distribution <- ((k/0.05)+100)
} else {
  Distribution <- 1000
}

presso<-MRPRESSO::mr_presso(NbDistribution = Distribution,BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data =Mydata,  SignifThreshold = 0.05)
presso
########################################################################################
if ((k / 10) > 6) {
  height <- (k / 10)
} else {
  height <- 6
}

single <- mr_leaveoneout(Mydata)
my_plot<-mr_leaveoneout_plot(single)   #留一法检验敏感性
my_plot
ggplot2::ggsave("mr_leaveoneout_plot.png", plot = my_plot[[1]], dpi = 200, width = 8, height = height,limitsize = FALSE)

my_plot<-mr_scatter_plot(res,Mydata)#散点图
ggplot2::ggsave("mr_scatter_plot.png", plot = my_plot[[1]], dpi = 200, width = 8, height = 6)

res_single <- mr_singlesnp(Mydata)
my_plot<-mr_forest_plot(res_single)#森林图
ggplot2::ggsave("mr_forest_plot.png", plot = my_plot[[1]], dpi = 200, width = 8, height = height,limitsize = FALSE)

my_plot<-mr_funnel_plot(res_single)#漏斗图
ggplot2::ggsave("mr_funnel_plot.png", plot = my_plot[[1]], dpi = 200, width = 8, height = 6)


# ###############################################################steiger分析
# #增加样本量数据
# n <- which(colnames(Mydata) == "samplesize.outcome")
# Mydata[n] <- 256523  # 将samplesize.outcome所在列所有值替换为132456
# m <- which(colnames(Mydata) == "samplesize.exposure")
# Mydata[m] <- 246139  # 将samplesize.exposure所在列所有值替换为132456
# Mydata$samplesize.exposure <- 246139
# Mydata$samplesize.outcome <- 135430
#############################
Mydata$r.exposure<-get_r_from_bsen(b = Mydata$beta.exposure,se = Mydata$se.exposure,n = Mydata$samplesize.exposure)
Mydata$r.outcome<-get_r_from_bsen(b = Mydata$beta.outcome,se = Mydata$se.outcome,n = Mydata$samplesize.outcome)
######
# #单个steiger分析？
# TwoSampleMR::steiger_filtering(dat = Mydata)
# #详细的steiger
# a<-TwoSampleMR::mr_steiger(p_exp = Mydata$pval.exposure,
#                         p_out = Mydata$pval.outcome,
#                         n_exp = Mydata$samplesize.exposure,
#                         n_out = Mydata$samplesize.outcome,
#                         r_exp = Mydata$r.exposure,
#                         r_out = Mydata$r.outcome
#                         )
#简单的steiger
steiger_directionality_test<-TwoSampleMR::directionality_test(dat = Mydata)
steiger_directionality_test



#############################################################################

result<-generate_odds_ratios(res)#算OR值
result

# res$lci95 <-res$b-1.96*res$se
# res$upi95 <-res$b+1.96*res$se
# res

write.csv(result, file="res.csv")
write.csv(heterogeneity, file="heterogeneity.csv")
#write.csv(single, file="single.csv")
write.csv(pleio, file="pleio.csv")
write.csv(presso[["Main MR results"]], file="presso.csv")
write.csv(steiger_directionality_test, file="steiger_directionality_test.csv")
#################################################################################
# Extract required parameters
EAF=Mydata$eaf.exposure
beta=Mydata$beta.exposure
N=Mydata$samplesize.exposure
SE=Mydata$se.exposure
k=nrow(Mydata)
# ###################################################################小R2计算方法  : R^2 = 2 * EAF * (1-EAF) * beta^2/((SE^2)*N)
# Mydata$Rsquare = (2*EAF*(1-EAF)*(beta^2))/((SE^2)*N)
# Mydata$Rsquare
# #####################较小的F计算方法
# Mydata$Fstatistic=((N-k-1)*Mydata$Rsquare)/(k*(1-Mydata$Rsquare))
# Mydata$Fstatistic
# ################################################################### 大R2计算方法 一
# Mydata$Rsquare = 2*EAF*(1-EAF)*(beta^2)
# ####################################
# sum(Mydata$Rsquare)
# sum(Mydata$Fstatistic)
################################################################### 大R2计算方法 二
# # R^2 = 2 * EAF * (1-EAF) * beta^2 / 
##         [(2 * EAF * (1-EAF) * beta^2)+
##           (2 * EAF * (1-EAF) * N * SE(beta)^2)]
Mydata$Rsquare=2*EAF*(1-EAF)*(beta^2) / (2*EAF*(1-EAF)*(beta^2)+2*EAF*(1-EAF)*N*(SE^2))

#####################较小的F计算方法
Mydata$Fstatistic=((N-k-1)*Mydata$Rsquare)/(k*(1-Mydata$Rsquare))
Mydata$Fstatistic
# # #####################较大的F计算方法
# Mydata$Fstatistic=((N-2)*Mydata$Rsquare)/(1-Mydata$Rsquare)
# Mydata$Fstatistic

sum(Mydata$Rsquare)
sum(Mydata$Fstatistic)
Mydata$sumR2<-sum(Mydata$Rsquare)
Mydata$sumF<-sum(Mydata$Fstatistic)
####################################
##############################################################
write.csv(Mydata, file="mydata1.csv")

