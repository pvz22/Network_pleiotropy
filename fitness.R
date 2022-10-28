#Comparing fitness estimates to pleiotropy estimates for deletion strains, and comparing fitness of trans-regulatory deletions
#to fitness of cis-regulatory deletions
rm(list = ls())
library(dplyr)
library(ggplot2)
library(cowplot)

figdir <- "/Users/petravandezande/Documents/Figures/Projects/Pleiotropy/Friend/Revision"
outputdir <- "/Users/petravandezande/Documents/Output/Projects/Pleiotropy/Friend"
fitdatdir <- "/Users/petravandezande/Documents/Data/Projects/Pleiodel"
data.dir <- "/Users/petravandezande/Documents/Data/Projects/Pleiodel"
#GGplot theme
THEMEMAIN <- function() {
  theme_bw() +
    theme(axis.text = element_text(size = 25), axis.title = element_text(size = 25), plot.margin = unit(c(2,1,1,1),"cm"), plot.title = element_text(size = 25, hjust = 0.5))
}

Zhangdat <- read.table(paste0(fitdatdir,"/supplementary_data_6.txt"), header = TRUE, stringsAsFactors = FALSE)
transsizes <- read.table(paste0(outputdir,"/transsizeslim.txt"), sep = "\t", stringsAsFactors = FALSE, header = 1)
KS1 <- read.table(paste0(data.dir,"/KemmerendataS1.txt"), sep = "\t", header = 1)
KS1$gene <- tolower(KS1$gene)

missingvec <- c()
for (i in 1:nrow(transsizes)) {
  if(transsizes[i,"Strain"] %in% c(KS1$gene)) {
    transsizes[i,"Transstrain"] <- KS1[KS1$gene == transsizes[i,"Strain"],"orf.name"]
  }
  else {
    missingvec <- c(missingvec,transsizes[i,"Strain"])
  }
  if(transsizes[i,"Focalgene"] %in% c(KS1$gene)) {
    transsizes[i,"Cisstrain"] <- KS1[KS1$gene == transsizes[i,"Focalgene"],"orf.name"]
  }
  else {
    missingvec <- c(missingvec,transsizes[i,"Focalgene"])
  }
}
length(unique(missingvec)) #There are 15 genes named with aliases that are different in this table

for (i in 1:nrow(transsizes)) {
  if(is.na(transsizes[i,"Transstrain"])) {
    transsizes[i,"Transstrain"] <- toupper(transsizes[i,"Strain"])
  }
  if(is.na(transsizes[i,"Cisstrain"])) {
    transsizes[i,"Cisstrain"] <- toupper(transsizes[i,"Focalgene"])
  }
}

transsizes[transsizes$Strain == "kem1","Transstrain"] <- "YGL173C"
transsizes[transsizes$Focalgene == "kem1","Cisstrain"] <- "YGL173C"
transsizes[transsizes$Focalgene == "lug1","Cisstrain"] <- "YLR352W"
transsizes[transsizes$Strain == "lug1","Transstrain"] <- "YLR352W"
transsizes$Cisstrain <- as.character(transsizes$Cisstrain)
transsizes$Transstrain <- as.character(transsizes$Transstrain)

for (i in 1:nrow(transsizes)) {
  try(
    transsizes[i,"Transfit"] <- Zhangdat[Zhangdat$Strain == transsizes[i,"Transstrain"],"Fitness_relative_to_HO_in_YPD"])
  try(
    transsizes[i,"Cisfit"] <- Zhangdat[Zhangdat$Strain == transsizes[i,"Cisstrain"],"Fitness_relative_to_HO_in_YPD"])
}

nofitcis <- c(unique(transsizes[is.na(transsizes$Cisfit),"Cisstrain"])) #double checking that these are indeed not in the fitness dataset
table(nofitcis %in% Zhangdat$Strain)
nofittrans <- c(unique(transsizes[is.na(transsizes$Transfit),"Transstrain"])) #Neither are these
table(nofittrans %in% Zhangdat$Strain)

ggplot(data = transsizes, aes(x = 1-Cisfit, y = 1-Transfit)) +
  geom_point(color = "darkorchid4") +
  THEMEMAIN() +
  geom_abline(slope = 1, intercept = 0, color = "seagreen4", size = 2) +
  xlab("Fitness cost of \ncis-regulatory mutation") +
  ylab("Fitness cost of \ntrans-regulatory mutations")
ggsave("Fitnessscatall.pdf", plot = last_plot(), path = figdir, width = 7, height = 7)

#Removing lines with NAs
transsizes <- transsizes[!is.na(transsizes$Cisfit),]
transsizes <- transsizes[!is.na(transsizes$Transfit),]

nrow(transsizes[transsizes$Transfit < transsizes$Cisfit,]) #3951/4985

trans_med_fit <- transsizes %>%
  group_by(Focalgene) %>%
  summarize(medfit = median(Transfit))
trans_med_fit <- as.data.frame(trans_med_fit)

for (i in 1:nrow(transsizes)) {
  transsizes[i,"trans_med_fit"] <- trans_med_fit[trans_med_fit$Focalgene == transsizes[i,"Focalgene"],"medfit"]
}

ggplot(data = transsizes, aes(x = (1-Cisfit), y = (1-trans_med_fit))) +
  geom_point(color = "black") +
  geom_abline(slope = 1, intercept = 0, color = "seagreen4", size = 2) +
  THEMEMAIN() +
  xlab("Fitness cost of\n cis-regulatory mutation") +
  ylab("Median fitness cost of\n trans-regulatory mutations")
ggsave("Fitnessscatmed.pdf", plot = last_plot(), path = figdir, width = 7, height = 7)

length(unique(transsizes[transsizes$trans_med_fit <= transsizes$Cisfit,"Focalgene"])) #489/583

#Overall relationship between fitness and pleiotropy
cissizes <- read.table(paste0(outputdir,"/cissizelim.txt"), sep = "\t", header = 1)
cissizes$gene <- rownames(cissizes)
KS1 <- left_join(KS1, cissizes, by = "gene")
KS1$Strain <- KS1$orf.name
KS1wfit <- left_join(KS1, Zhangdat, by = "Strain")

ggplot(data = KS1wfit[!is.na(KS1wfit$cissize),], aes(x = cissize, y = Fitness_relative_to_HO_in_YPD)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue2") +
  THEMEMAIN() +
  xlab("Pleiotropy\n(Number of DEGs)") +
  ylab("Fitness")
ggsave("Fitmod.pdf", plot = last_plot(), path = figdir, width = 7, height = 7)

fitmod <- lm(KS1wfit$Fitness_relative_to_HO_in_YPD ~ KS1wfit$SIGCHG)

#Now looing at the with the slow growth signature removed
Genes <- data.frame(Strain = transsizes$Strain, Eudist = transsizes$Eudisttranssub, Fitness = transsizes$Transfit)
FCGenes <- data.frame(Strain = transsizes$Focalgene, Eudist = transsizes$Eudistcissub, Fitness = transsizes$Cisfit)
Allgenes <- rbind(Genes, FCGenes)
Allgenes <- unique(Allgenes)

ggplot(data = Allgenes, aes(x = Eudist, y = Fitness)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm") +
  THEMEMAIN() +
  xlab("Pleiotropy\n (Euclidean distance)\n Slow growth removed") +
  ylab("Fitness")
ggsave("Fitmodeudsub.pdf", plot = last_plot(), path = figdir, width = 7, height = 7)

fitmodsub <- lm(Allgenes$Fitness ~ Allgenes$Eudist)


#Looking at the fitness costs of cis and trans-regulatory mutations for different network cutoffs
#Reading in the trans-regulatory info talbes for each 'paramsweep' network
sets <- c(11,12,13,14,15,21,22,23,24,25,31,32,33,34,35,41,42,43,44,45)
KS1wfit[KS1wfit$gene == "xrn1","gene"] <- "kem1"
KS1wfit[KS1wfit$gene == "ylr352w","gene"] <- "lug1"
KS1wfit[KS1wfit$gene == "coq8","gene"] <- "abc1"
KS1wfit[KS1wfit$gene == "yjr084w","gene"] <- "csn12"

for (i in 1:length(sets)) {
  x <- read.table(paste0(outputdir,"/transsizes",sets[i],".txt"), sep = "\t", header = 1)
  assign(paste0("transsizes",sets[i]),x)
}

for (n in 1:length(sets)) {
  transtemp <- get(paste0("transsizes",sets[n]))
  for (i in 1:nrow(transtemp)) {
    if (transtemp[i,"Focalgene"] %in% KS1wfit$gene) {
      transtemp[i,"Cisfit"] <- KS1wfit[KS1wfit$gene == transtemp[i,"Focalgene"],"Fitness_relative_to_HO_in_YPD"]
    }
    else {
      transtemp[i,"Cisfit"] <- KS1wfit[tolower(KS1wfit$orf.name) == transtemp[i,"Focalgene"],"Fitness_relative_to_HO_in_YPD"]
    }
    if (transtemp[i,"Strain"] %in% KS1wfit$gene) {
      transtemp[i,"Transfit"] <- KS1wfit[KS1wfit$gene == transtemp[i,"Strain"],"Fitness_relative_to_HO_in_YPD"]
    }
    else (
      transtemp[i,"Transfit"] <- KS1wfit[tolower(KS1wfit$orf.name) == transtemp[i,"Strain"],"Fitness_relative_to_HO_in_YPD"]
    )
  }
  assign(paste0("transsizes",sets[n]), transtemp)
}

for (i in 1:length(sets)) {
  transtemp <- get(paste0("transsizes",sets[i]))
  transtemp <- data.frame(transtemp, stringsAsFactors = FALSE)
  trans_med <- transtemp %>%
    group_by(Focalgene) %>%
    summarize(med = median(Transfit, na.rm = TRUE))
  trans_med <- as.data.frame(trans_med)
  for (j in 1:nrow(transtemp)) {
    transtemp[j,"Medfit"] <- trans_med[as.character(trans_med$Focalgene) == as.character(transtemp[j,"Focalgene"]),"med"]
  }
  plot <- ggplot(data = transtemp, aes(x=1-Cisfit, y = 1-Medfit)) +
    geom_point(color = "black") +
    theme_bw() +
    geom_abline(slope = 1, intercept = 0, color = "seagreen4", size = 2) +
    xlab("") +
    ylab("")
  assign(paste0("p",sets[i]), plot)
}
plotlist <- list(p11,p12,p13,p14,p15,p21,p22,p23,p24,p25,p31,p32,p33,p34,p35,p41,p42,p43,p44,p45)
plot_grid(plotlist = plotlist, nrow = 4, ncol = 5)
ggsave("paramsweepfit.png",plot = last_plot(), path = figdir, width = 12, height = 12)

#Now looking at fitness in other environments
Zhangdatetoh <- read.table(paste0(fitdatdir,"/maclean_fit_Etoh.txt"), header = FALSE, skip = 1, stringsAsFactors = FALSE)
Zhangdath2o2 <- read.table(paste0(fitdatdir,"/maclean_fit_H2O2.txt"), header = FALSE, skip = 1, stringsAsFactors = FALSE)
Zhangdatnacl <- read.table(paste0(fitdatdir,"/maclean_fit_NaCl.txt"), header = FALSE, skip = 1, stringsAsFactors = FALSE)
Zhangdatcocl2 <- read.table(paste0(fitdatdir,"/maclean_fit_CoCl2.txt"), header = FALSE, skip = 1, stringsAsFactors = FALSE)
Zhangdatso <- read.table(paste0(fitdatdir,"/maclean_fit_SO.txt"), header = FALSE, skip = 1, stringsAsFactors = FALSE)
Zhangdattemp <- read.table(paste0(fitdatdir,"/Maclean_fit_temp.txt"), header = FALSE, skip = 1, stringsAsFactors = FALSE)

Zhangdattog <- cbind(Zhangdattemp[,1:3],Zhangdatetoh[,3],Zhangdath2o2[,3],Zhangdatnacl[,3],Zhangdatcocl2[,3],Zhangdatso[,3])
colnames(Zhangdattog) <- c("orf.name","YPD","Temp","Eth","Perox","NaCl","CoCl2","SO")

head(KS1)

TOG <- left_join(Zhangdattog, KS1, by = "orf.name")

a <- ggplot(data = TOG, aes(x = SIGCHG, y = Temp)) +
  geom_point() +
  geom_smooth(method = "lm") +
  THEMEMAIN() +
  xlab("Number of DEGs in SC") +
  ylab("Fitness at 40C")

b <- ggplot(data = TOG, aes(x = SIGCHG, y = Eth)) +
  geom_point() +
  geom_smooth(method = "lm") +
  THEMEMAIN() +
  xlab("Number of DEGs in SC") +
  ylab("Fitness in Ethanol")

c <- ggplot(data = TOG, aes(x = SIGCHG, y = Perox)) +
  geom_point() +
  geom_smooth(method = "lm") +
  THEMEMAIN() +
  xlab("Number of DEGs in SC") +
  ylab("Fitness in Hydrogen Peroxide")

d <- ggplot(data = TOG, aes(x = SIGCHG, y = NaCl)) +
  geom_point() +
  geom_smooth(method = "lm") +
  THEMEMAIN() +
  xlab("Number of DEGs in SC") +
  ylab("Fitness in high NaCl")

e <- ggplot(data = TOG, aes(x = SIGCHG, y = CoCl2)) +
  geom_point() +
  geom_smooth(method = "lm") +
  THEMEMAIN() +
  xlab("Number of DEGs in SC") +
  ylab("Fitness in Cobalt chloride")

f <- ggplot(data = TOG, aes(x = SIGCHG, y = SO)) +
  geom_point() +
  geom_smooth(method = "lm") +
  THEMEMAIN() +
  xlab("Number of DEGs in SC") +
  ylab("Fitness in Superoxide")

plot_grid(a,b,c,d,e,f, nrow = 2)
ggsave("Fitnessenvs.png", plot = last_plot(), path = figdir, width = 20, height = 15)
