# Rscript for analysis of balsam poplar common garden height increment growth as a function of genetic offsets and climate distance
# Accompanies: Fitzpatrick, M.C., V.E. Chhatre, R. Soolanayakanahally, and S.R. Keller (in press). Experimental support for genomic prediction of climate maladaptation using the machine learning approach Gradient Forests. Molecular Ecology Resources. 
# Authorea preprint: https://doi.org/10.22541/au.159863198.86187354
# Script by: Stephen R. Keller (srkeller@uvm.edu); Dept. of Plant Biology; University of Vermont
# Last modified: Mar 02, 2021

library(lmerTest)
library(ggplot2)
library(sjPlot)
library(MuMIn)
library(cowplot)
library(grid)
library(gridExtra)

# Outlier SNPs:
GF_predictSTD <- read.csv("core336.gf.stdFreqs.gardenOffset.csv") # These are based on using the XTX standardized allele freqs before feeding to GF
colnames(GF_predictSTD) <- c("pop_code", "VT", "IH")
GF_predictAllOutliers <- read.csv("core336.gf.allOls.gardenOffset.csv") # These are based on the concatenated set of all outlier loci from Bayenv, LFMM, and native GF analyses
colnames(GF_predictAllOutliers) <- c("pop_code", "VT", "IH")

# Random SNPs:
GF_predictRAND <- read.csv("core336.500_RANDOM.gardenOffset.csv")
colnames(GF_predictAllOutliers) <- c("pop_code", "VT", "IH")

GF_predictRANDavg <- aggregate(GF_predictRAND, by=list(GF_predictRAND$pop_code), FUN="mean")

Clim <- read.csv("core336.mahal.gardenOffset.csv")

pairs(cbind(GF_predictAllOutliers[,2:3],GF_predictSTD[,2:3],GF_predictRANDavg[,3:4],Clim[,2:3]), 
        labels = c("All outliers: VT","All outliers: IH","GF-X: VT","GF-X: IH","Random: VT","Random: IH","Climate-only: VT","Climate-only: IH"),
        col="red",pch=19)

cor.test(GF_predictSTD[,2],GF_predictSTD[,3]) # corr of GF-X across gardens
cor.test(Clim[,2],Clim[,3]) # corr of clim-only across gardens

# Poplar common garden height data:
load("gardenQuery.RData")

HI_VT_IH <- merge(HI_VT, joiner)

# Remove missing obs
HI_VT_IH[HI_VT_IH$height_increment==0,4] <- NA # Set missing height data from '0' to 'NA'
length(HI_VT_IH[which(is.na(HI_VT_IH$height_increment)==F),4])  #1461 obs

# Join with GF genetic offset predictions
GF_HI_VT_IH <- merge(GF_predictSTD[,1:2], HI_VT_IH, by.x="pop_code", by.y="pop_code")
names(GF_HI_VT_IH)[2] = "Offset"
GF_HI_VT_IH$Offset2 = GF_HI_VT_IH$Offset^2

# Use MLM to factor out microsite effects within each garden and then model pop_code BLUPs as a function of Offset

VT_MLM <- lmer(height_increment ~ tree_num*row_num + (1|pop_code) + (1|ind_code), data=subset(HI_VT_IH, garden_id=="VT"))
summary(VT_MLM)

VTblups <- coef(VT_MLM)[2] # Gets vector of BLUPs for the populations
pop_code <- row.names(VTblups$pop_code)

VThtblups <- data.frame(pop_code, VTblups[1])
colnames(VThtblups) = c("pop_code","HtBlups")
VThtblups$Garden=as.factor("VT")

# Merge height BLUPs with GF genetic offsets based on GF-X, All outliers, Random SNPs, and climate-only (respectively)
GFX_VThtblups <- merge(GF_predictSTD[,c(1,2)], VThtblups[,c(1,2,6)], by="pop_code") 
GFAll_VThtblups <- merge(GF_predictAllOutliers[,c(1,2)], VThtblups[,c(1,2,6)], by="pop_code")
GFRand_VThtblups <- merge(GF_predictRANDavg[,c(1,3)], VThtblups[,c(1,2,6)], by.x="Group.1",by.y="pop_code")
Clim_VThtblups <- merge(Clim[,c(1,2)], VThtblups[,c(1,2,6)], by="pop_code")

colnames(GFX_VThtblups) = c("pop_code","offset", "HtBlups", "Garden")
colnames(GFAll_VThtblups) = c("pop_code","offset", "HtBlups", "Garden")
colnames(GFRand_VThtblups) = c("pop_code","offset", "HtBlups", "Garden")
colnames(Clim_VThtblups) = c("pop_code","offset", "HtBlups", "Garden")


HI_IH <- read.csv("IH_phenology_2015.csv", header=T)

IH_MLM2 <- lmer(height_increment ~ field + (1|pop_code) + (1|ind_code), data=HI_IH[which(HI_IH$field != "P1" & HI_IH$field != "S2"),])
summary(IH_MLM2)

IHblups <- coef(IH_MLM2)[2] # Gets vector of BLUPs for the populations
pop_code <- row.names(IHblups$pop_code)

IHhtblups <- data.frame(pop_code, IHblups)
colnames(IHhtblups) = c("pop_code","HtBlups")
IHhtblups$Garden=as.factor("IH")

GFX_IHhtblups <- merge(GF_predictSTD[,c(1,3)], IHhtblups[,c(1,2,4)], by="pop_code")
GFAll_IHhtblups <- merge(GF_predictAllOutliers[,c(1,2)], IHhtblups[,c(1,2,4)], by="pop_code")
GFRand_IHhtblups <- merge(GF_predictRANDavg[,c(1,4)], IHhtblups[,c(1,2,4)], by.x="Group.1",by.y="pop_code")
Clim_IHhtblups <- merge(Clim[,c(1,3)], IHhtblups[,c(1,2,4)], by="pop_code")

colnames(GFX_IHhtblups) = c("pop_code","offset", "HtBlups", "Garden")
colnames(GFAll_IHhtblups) = c("pop_code","offset", "HtBlups", "Garden")
colnames(GFRand_IHhtblups) = c("pop_code","offset", "HtBlups", "Garden")
colnames(Clim_IHhtblups) = c("pop_code","offset", "HtBlups", "Garden")

# Merge height BLUPs and GF genetic offset values across both common gardens
GFX_Combined <- rbind(GFX_VThtblups,GFX_IHhtblups)
GFAll_Combined <- rbind(GFAll_VThtblups,GFAll_IHhtblups)
GFRand_Combined <- rbind(GFRand_VThtblups,GFRand_IHhtblups)
Clim_Combined <- rbind(Clim_VThtblups,Clim_IHhtblups)


GFX_mod <- lm(HtBlups ~ offset + I(offset^2) + Garden, data=GFX_Combined)
summary(GFX_mod)

GFX_plot <- plot_model(GFX_mod, type="pred", terms=c("offset"), show.data=T, axis.title=c("Genetic offset", "Height increment (cm)"), title="",dot.size = 4,font_size(4))
GFX_plot

GFAll_mod <- lm(HtBlups ~ offset + I(offset^2) + Garden, data=GFAll_Combined)
summary(GFAll_mod)

GFAll_plot <- plot_model(GFAll_mod, type="pred", terms=c("offset"), show.data=T, axis.title=c("Genetic offset", "Height increment (cm)"), title="",dot.size = 4,font_size(4))
GFAll_plot

GFRand_mod <- lm(HtBlups ~ offset + I(offset^2) + Garden, data=GFRand_Combined)
summary(GFRand_mod)

GFRand_plot <- plot_model(GFRand_mod, type="pred", terms=c("offset"), show.data=T, axis.title=c("Genetic offset", "Height increment (cm)"), title="",dot.size = 4,font_size(4))
GFRand_plot

Clim_Combined_mod <- lm(HtBlups ~ offset + Garden, data=Clim_Combined)
summary(Clim_Combined_mod)

Clim_Combined_plot <- plot_model(Clim_Combined_mod, type="pred", terms=c("offset"), show.data=T, axis.title=c("Climate-only offset", "Height increment (cm)"), title="",dot.size = 4,font_size(4))
Clim_Combined_plot

# Calculate AIC scores for model comparison
AICc(GFX_mod) #GF-X
AICc(GFAll_mod) #All
AICc(GFRand_mod) #Rand
AICc(Clim_Combined_mod) #CLim

# Compare GF-R2 for SNPs that are unique versus shared outliers among the different sel'n scan methods:
R2 <- read.csv("snpR2.csv", header=T)

freqR2 <- rowSums(table(R2$SNP,R2$model))

SNP <- names(freqR2)

df <- as.data.frame(cbind(SNP,freqR2))

R2_v2 <- merge(R2, df, by="SNP")

p5 <- ggplot(R2_v2, aes(x=model, y=R2,fill=model)) + 
  geom_violin(trim=FALSE) +
  scale_x_discrete(limits=c("GF-Raw", "GF-X", "LFMM", "Bayenv")) +
  labs(x = "Outlier Test", y = expression(paste("GF Model ",R^2))) +
  theme(legend.position = "none")

p5 + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.25)

modelR2 <- lm(R2~model, data=R2_v2)
aovR2 <- aov(modelR2)
TukeyHSD(aovR2)

p6 <- ggplot(R2_v2, aes(x=freqR2, y=R2,fill=freqR2)) + 
  geom_violin(trim=FALSE) +
  labs(x="Number of overlapping tests", 
       y = expression(paste("GF Model ",R^2))) + 
       theme(legend.position = "none")

p6 + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.35)


testsR2 <- lm(R2~as.factor(freqR2), data=R2_v2)

aovtestsR2 <- aov(testsR2)
TukeyHSD(aovtestsR2)

