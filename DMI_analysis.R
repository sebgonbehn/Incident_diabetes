setwd("/Users/sebastianbehn/Downloads/")
library(plyr);library(readr);library(readxl);library(ggplot2);library(data.table);library(magrittr);library(dplyr);library(ggpubr);library(corrplot);library(reshape2);
library(car);library(patchwork);library(MASS);library(ggridges);library(ggstance);library(haven);library(survival);library(survminer);library(survivalROC);
library(performance);library(rms);library(caret);library(mice);library(VIM);library(glmnet);library(timeROC);library(ggridges);library(forcats)

#Dataset loading
base <- read_sav(file.path("cohorte_final.sav"))
base <- base[base$pac_loc=="1",]
ids1 <- read_xls("ids_UIEM-UIDS.xls")
colnames(ids1)[which(names(ids1) == "id_endocrin")] <- "id_uiem"
base <- merge(ids1,base,by="id_uiem",all.y=T)
ipac <- read_xlsx("ipaqbasal.xlsx")
base <- merge(base,ipac,by="id_uiem")

#Add genetic variables
geno <- read.table(file="all_pheno.txt",he=T)
ids <- read_xlsx("IDsEquivalentsUIEM-Paivi.xlsx")
metsb <- merge(geno,ids,by="IID")
metsb <- merge(metsb,base,by="id_uiem")
metsb$iid <- NULL
colnames(base)[which(names(base) == "iid")] <- "IID"
sigma <- merge(geno,base,by="IID")
sigma <- sigma[!sigma$id_uiem%in%metsb$id_uiem,]
data <- rbind(metsb,sigma)
prs <- read.table(file="sigma_metsb_prs.txt",he=T)
prs$Cohort <- NULL
pcs <- read.table(file="PCs_GA_aims_sigma_metsb_covid.txt",he=T)
colnames(pcs)[which(names(pcs) == "AMR")] <- "AMR_aims"
colnames(pcs)[which(names(pcs) == "EUR")] <- "EUR_aims"
ga_pcs_metsb <- read.table(file="GlobalAnc/Metsb_PCs+GlobalAnc_bothGW.txt",he=T)
ga_pcs_metsb$group <- "METSB"
ga_pcs_omni <- read.table(file="GlobalAnc/OmniExchip_PCs+GlobalAnc_bothGW.txt",he=T)
ga_pcs_omni$FID.x <- NULL
ga_pcs_omni$group <- "OMNI"
ga_pcs_refpanel <- read.table(file="GlobalAnc/RefPanel4k_PCs+GlobalAnc_bothGW.txt",he=T)
ga_pcs_refpanel$FID.x <- NULL
ga_pcs_refpanel$group <- "REFP"
ga_pcs_sigma3 <- read.table(file="GlobalAnc/Sigma3_PCs+GlobalAnc_bothGW.txt",he=T)
colnames(ga_pcs_sigma3)[which(names(ga_pcs_sigma3) == "ID_1")] <- "IID"
ga_pcs_sigma3$group <- "SIGMA3"
ga_pcs <- rbind(ga_pcs_metsb,ga_pcs_omni,ga_pcs_refpanel,ga_pcs_sigma3)
pcs <- merge(pcs,ga_pcs,by="IID")
data <- merge(prs,data,by="IID")
data <- merge(pcs,data,by="IID")

#Data preparation
data$age2 <- data$age^2
data$tiempo_anios <- as.numeric(as.Date(data$fechamuestra2,"%Y-%m-%d") - as.Date(data$fechamuestra1,"%Y-%m-%d"))/365
data$tiempo_anios <- ifelse(data$tiempo_anios<0,NA,data$tiempo_anios)
data <- within(data, {
  PRS_T2D[scale(PRS_T2D)> 3] <- NA
  PRS_T2D[scale(PRS_T2D)< -3] <- NA
  PC1_aims[scale(PRS_T2D)> 3] <- NA
  PC1_aims[scale(PRS_T2D)< -3] <- NA
  años_est}) 
data <- data[!is.na(data$PRS_T2D),]
data <- data[!is.na(data$PC1_aims),]
data <- data[!is.na(data$tiempo_anios),]
data <- data[!is.na(data$dm_incidente),]
data$PRS_T2Dadj <- summary(lm(PRS_T2D~AMR,data=data))$res
data$PRS_BMIadj <- summary(lm(PRS_BMI~AMR,data=data))$res
data$Cohort <- ifelse(data$Cohort!="METSB","SIGMA","METSB")

#Variable definition
data <- within(data,{
  obe_inc <- NA
  obe_inc[imc_basal<30&imc_final<30] <- 0
  obe_inc[imc_basal<30&imc_final>=30] <- 1})
data$icc1 <- data$cintura1 / data$cadera1
data$PRS_T2Dadjstd <- scale(data$PRS_T2Dadj)
data$PRS_BMIadjstd <- scale(data$PRS_BMIadj)

#Generate density and scatter plots
a <- ggplot(data[!is.na(data$dm_incidente),],aes(x = AMR*100, y = as.factor(dm_incidente),fill=as.factor(dm_incidente))) +labs(y="Incident T2D",x="NAT Ancestry (%)")+scale_fill_manual(values=c("salmon4","salmon"))+
  stat_density_ridges(quantile_lines = TRUE, alpha = 0.5,bandwidth = 3.5,scale = 1.2)+scale_y_discrete(labels=c("Controls","Cases"))+theme_minimal()+guides(fill=F)

b <- ggplot(data[!is.na(data$dm_incidente),],aes(x=as.factor(dm_incidente),y=scale(PRS_T2Dadj)))+geom_boxplot(alpha=0.5,aes(fill=as.factor(dm_incidente),group=as.factor(dm_incidente)),width=0.3,size = 0.2)+geom_violin(alpha=0.3,aes(fill=as.factor(dm_incidente),group=as.factor(dm_incidente)),size = 0.05,adjust=0.4) + 
  ylab("PRS T2D")+xlab("Incident T2D")+scale_x_discrete(labels=c("Controls","Cases"))+
  stat_summary(fun="mean", geom="point", shape=23, size=2, fill="grey",alpha=0.8) +
  theme_minimal()+guides(fill=F)+scale_fill_manual(values=c("navajowhite2","navajowhite4"))

b <- ggplot(data[!is.na(data$dm_incidente),],aes(x = scale(PRS_T2Dadj), y = as.factor(dm_incidente),fill=as.factor(dm_incidente))) +labs(y="Incident T2D",x=" PRS T2D")+scale_fill_manual(values=c("navajowhite4","navajowhite"))+
  stat_density_ridges(quantile_lines = TRUE, alpha = 0.5,bandwidth = 0.18,scale = 1.2)+scale_y_discrete(labels=c("Controls","Cases"))+theme_minimal()+guides(fill=F)

c <- ggplot(data[!is.na(data$glucosa1),], aes(x=AMR*100, y=glucosa1)) + 
  geom_point(shape=18, color="gray80",size=1.5)+labs(y="Basal glucose (mg/dL)",x="NAT Ancestry (%)")+
  geom_smooth(method=lm,color="salmon2", fill="salmon2")+theme_minimal()+ stat_cor(method = "pearson", label.y = 130, label.x = 0)
d <- ggplot(data[!is.na(data$glucosa1),], aes(y=glucosa1, x=scale(PRS_T2Dadj))) + 
  geom_point(shape=18, color="gray80",size=1.5)+labs(y="Basal glucose (mg/dL)",x="PRS T2D")+
  geom_smooth(method=lm,color="navajowhite2", fill="navajowhite2")+theme_minimal()+ stat_cor(method = "pearson", label.y = 130, label.x = -3)

e <- a/c+plot_layout(heights = c(1,2)) & scale_x_continuous(limits=c(0,100))
f <- b/d+plot_layout(heights = c(1,2)) & scale_x_continuous(limits=c(-3,3))
plot <- e|f

a <- ggplot(data[!is.na(data$obe_inc),],aes(x = AMR*100, y = as.factor(obe_inc),fill=as.factor(obe_inc))) +labs(y="Incident Obesity",x="NAT Ancestry (%)")+scale_fill_manual(values=c("salmon4","salmon"))+
  stat_density_ridges(quantile_lines = TRUE, alpha = 0.5,bandwidth = 3.5,scale = 1.2)+scale_y_discrete(labels=c("Controls","Cases"))+theme_minimal()+guides(fill=F)

b <- ggplot(data[!is.na(data$obe_inc),],aes(x = scale(PRS_BMIadj), y = as.factor(obe_inc),fill=as.factor(obe_inc))) +labs(y="Incident Obesity",x=" PRS BMI")+scale_fill_manual(values=c("navajowhite4","navajowhite"))+
  stat_density_ridges(quantile_lines = TRUE, alpha = 0.5,bandwidth = 0.18,scale = 1.2)+scale_y_discrete(labels=c("Controls","Cases"))+theme_minimal()+guides(fill=F)

c <- ggplot(data[!is.na(data$imc_basal),], aes(x=AMR*100, y=imc_basal)) + ylim(15,55)+
  geom_point(shape=18, color="gray80",size=1.5)+labs(y=bquote('BMI '~(kg/m^2)),x="NAT Ancestry (%)")+
  geom_smooth(method=lm,color="salmon2", fill="salmon2")+theme_minimal()+ stat_cor(method = "pearson", label.y = 50, label.x = 0)
d <- ggplot(data[!is.na(data$imc_basal),], aes(y=imc_basal, x=scale(PRS_BMIadj))) + 
  geom_point(shape=18, color="gray80",size=1.5)+labs(y=bquote('BMI '~(kg/m^2)),x="PRS BMI")+ylim(15,55)+
  geom_smooth(method=lm,color="navajowhite2", fill="navajowhite2")+theme_minimal()+ stat_cor(method = "pearson", label.y = 50, label.x = -3)

e <- a/c+plot_layout(heights = c(1,2)) & scale_x_continuous(limits=c(0,100))
f <- b/d+plot_layout(heights = c(1,2)) & scale_x_continuous(limits=c(-3,3))
tiff("test.tiff", units="in", width=8, height=7, res=1200)
plot2 <- e|f
plot/plot2
dev.off()

####Cox models

##Incident T2D
effect <- NA
eff <- as.data.frame(summary(coxph(Surv(tiempo_anios, dm_incidente) ~ scale(edad)+as.factor(sexo)+scale(AMR)+scale(PRS_T2Dadjstd), data = data))$coefficients)
eff$Variable <- factor(c("Age","Sex (male)","NAT","PRS T2D"))
eff$LCL <- eff$`exp(coef)` - 1.96*eff$`se(coef)`
eff$UCL <- eff$`exp(coef)` + 1.96*eff$`se(coef)`
eff$Sign <- ifelse(eff$`Pr(>|z|)` <= 0.01,1,0)
eff$P <- signif(eff$`Pr(>|z|)`,digits=2)
eff$P <- paste0("P = ",eff$P)
eff$Response <- "Incident T2D"
effect <- rbind(eff,effect)
eff <- na.omit(effect)

forest1 <- eff %>%mutate(Variable = factor(Variable, levels=c("PRS T2D", "NAT","Sex (male)", "Age")))%>%
  ggplot(., aes(y = exp(coef), x = Variable, color=factor(Sign))) +
  stat_summary(fun="mean", geom="errorbar") + 
  geom_pointrange(aes(ymax = UCL, ymin = LCL))+
  geom_point(size = 3,shape=18)+coord_flip()+ylim(c(0.45,1.45))+
  scale_color_manual(values=c("azure4", "brown2"))+
  ylab("Hazard ratio (95% CI)") +xlab("")+ggtitle("Incident T2D")+
  geom_hline(yintercept = 1, linetype = 3)+
  theme_bw()+theme(legend.position="none")+geom_text(aes(label=P),hjust=0.4, vjust=-0.8,size=3)+ 
  theme(plot.title = element_text(size = 13, face = "bold"))
tiff("test.tiff", units="in", width=7, height=4, res=1000)
forest1
dev.off()

##Incident obesity
effect <- NA
eff <- as.data.frame(summary(coxph(Surv(tiempo_anios, obe_inc) ~ scale(edad)+as.factor(sexo)+scale(AMR)+scale(PRS_BMIadj), data = data))$coefficients)
eff$Variable <- factor(c("Age","Sex (male)","NAT","PRS BMI"))
eff$LCL <- eff$`exp(coef)` - 1.96*eff$`se(coef)`
eff$UCL <- eff$`exp(coef)` + 1.96*eff$`se(coef)`
eff$Sign <- ifelse(eff$`Pr(>|z|)` <= 0.05,1,0)
eff$P <- signif(eff$`Pr(>|z|)`,digits=2)
eff$P <- paste0("P = ",eff$P)
eff$Response <- "Incident BMI"
effect <- rbind(eff,effect)
eff <- na.omit(effect)

forest2 <- eff %>%mutate(Variable = factor(Variable, levels=c("PRS BMI", "NAT","Sex (male)", "Age")))%>%
  ggplot(., aes(y = exp(coef), x = Variable, color=factor(Sign))) +
  stat_summary(fun="mean", geom="errorbar") + 
  geom_pointrange(aes(ymax = UCL, ymin = LCL))+
  geom_point(size = 3,shape=18)+coord_flip()+
  scale_color_manual(values=c("azure4", "brown2"))+
  ylab("Hazard ratio (95% CI)") +xlab("")+ggtitle("Incident obesity")+
  geom_hline(yintercept = 1, linetype = 3)+ylim(c(0.45,1.45))+
  theme_bw()+theme(legend.position="none")+geom_text(aes(label=P),hjust=0.4, vjust=-0.8,size=3)+ 
  theme(plot.title = element_text(size = 13, face = "bold"))

tiff("test.tiff", units="in", width=6, height=5, res=1200)
forest1/forest2
dev.off()

#Survival curves
data$quant_prs<- ifelse(data$PRS_T2Dadjstd>=mean(data$PRS_T2Dadjstd),1,0)
data$quant_amr<- ifelse(data$AMR>=mean(data$AMR),1,0)

tiff("test.tiff", units="in", width=6, height=9, res=1200)
a <- ggsurvplot(fit=survfit(Surv(tiempo_anios, dm_incidente)~quant_prs, data = data),data=data,title="T2D PRS stratification",pval = T, conf.int = T,
                risk.table = T,ggtheme = theme_minimal(),ncensor.plot=T,palette=c("salmon2","navajowhite4"),
                legend.labs = c("Low T2D PRS", "High T2D PRS"),ylim=c(0.4,1),xlim=c(0,5))
b <- ggsurvplot(fit=survfit(Surv(tiempo_anios, dm_incidente)~quant_amr, data = data),data=data,title="NAT ancestry stratification",pval = T, conf.int = T,
                risk.table = T,ggtheme = theme_minimal(),ncensor.plot=T,palette =c("salmon2","navajowhite4"),
                legend.labs = c("Low NAT ancestry", "High NAT ancestry"),ylim=c(0.4,1),xlim=c(0,5))
a
b
dev.off()



