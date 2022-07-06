rm(list = ls())
library(dplyr)
library(ggplot2)

figdir <- "/Users/petravandezande/Documents/Figures/Projects/Pleiotropy/Friend"
outputdir <- "/Users/petravandezande/Documents/Output/Projects/Pleiotropy/Friend"
fitdatdir <- "/Users/petravandezande/Documents/Data/Projects/Pleiodel"
THEMEMAIN <- function() {
  theme_bw() +
    theme(axis.text = element_text(size = 25), axis.title = element_text(size = 25), plot.margin = unit(c(2,1,1,1),"cm"), plot.title = element_text(size = 25, hjust = 0.5))
}

datacheck <- read.table(paste0(fitdatdir,"/dbxref.txt"), fill = NA, stringsAsFactors = FALSE)
Zhangdat <- read.table(paste0(fitdatdir,"/supplementary_data_6.txt"), header = TRUE, stringsAsFactors = FALSE)
transsizes <- read.table(paste0(outputdir,"/transsizeslim.txt"), sep = "\t", stringsAsFactors = FALSE)

alltrans <- data.frame(Strain = unique(as.character(transsizes$Strain)), stringsAsFactors = FALSE)
missingvec <- c()
for (i in 1:nrow(alltrans)) {
  if (!alltrans[i,"Strain"] %in% c(tolower(datacheck$V7))) {
    missingvec <- c(missingvec, alltrans[i,"Strain"])
  }
  else {
    try(
      alltrans[i,"FitnessJZstrain"] <- Zhangdat[Zhangdat$Strain == (datacheck[tolower(datacheck$V7) == alltrans[i,"Strain"],"V5"][1]),"Fitness_relative_to_HO_in_YPD"]
    )
  }
}
missingdf <- data.frame(Oldname = missingvec, stringsAsFactors = FALSE)
missingdf$Strain <- c("mot2","srb8","NA","ssn3","NA","xrn1",
                      "nut1","pgd1","srb5","elp3","gal11","sin4",
                      "cse2","ssn2","soh1","srb2","NA",'NA',"NA",
                      "iba57","gep3","rlf2","NA","bol3","gpb1","bol2","snx4",
                      "gpb2","dma1","coq8")
for (i in 1:nrow(alltrans)) {
  if(!alltrans[i,"Strain"] %in% c(tolower(datacheck$V7))) {
    x <- missingdf[missingdf$Oldname == alltrans[i,"Strain"],"Strain"]
    try(
      alltrans[i,"FitnessJZstrain"] <- Zhangdat[Zhangdat$Strain == (datacheck[tolower(datacheck$V7) == x,"V5"][1]),"Fitness_relative_to_HO_in_YPD"]
    )
  }
}

transsizes <- left_join(transsizes, alltrans, by = "Strain")
transsizes <- transsizes[transsizes$Strain != transsizes$Focalgene,]

trans_med <- transsizes %>%
  group_by(Focalgene) %>%
  summarize(med = median(FitnessJZstrain, na.rm = TRUE))

cissize <- read.table(paste0(outputdir,"/cissizelim.txt"), sep = "\t", stringsAsFactors = FALSE)
for (i in 1:nrow(cissize)) {
  try(
    cissize[i,"FitnessJZstrain"] <- Zhangdat[Zhangdat$Strain == (datacheck[tolower(datacheck$V7) == rownames(cissize)[i],"V5"][1]),"Fitness_relative_to_HO_in_YPD"]
  )
} 
for (i in c(which(rownames(cissize) %in% missingdf$Oldname))) {
  x <- missingdf[missingdf$Oldname == rownames(cissize)[i],"Strain"]
  try(
    cissize[i,"FitnessJZstrain"] <- Zhangdat[Zhangdat$Strain == (datacheck[tolower(datacheck$V7) == x,"V5"][1]),"Fitness_relative_to_HO_in_YPD"]
  )
}

for (i in 1:nrow(transsizes)) {
  transsizes[i,"cisfit"] <- cissize[rownames(cissize) == transsizes[i,"Focalgene"],"FitnessJZstrain"]
  transsizes[i,"trans_med"] <- trans_med[trans_med$Focalgene == transsizes[i,"Focalgene"],"med"]
}

ggplot(data = transsizes, aes(x = (1-cisfit), y = (1-trans_med))) +
  geom_point(color = "black") +
  geom_abline(slope = 1, intercept = 0, color = "seagreen4", size = 2) +
  THEMEMAIN() +
  xlab("Cis fitness cost\n(1-Fitness)") +
  ylab("Median trans fitness cost\n(1-Fitness)")
ggsave("Fitnessscatmed.pdf", plot = last_plot(), path = figdir, width = 7, height = 7)

ggplot(data = transsizes, aes(x = (1-cisfit), y = (1-FitnessJZstrain))) +
  geom_point(color = "black") +
  geom_abline(slope = 1, intercept = 0, color = "seagreen4", size = 2) +
  THEMEMAIN() +
  xlab("Cis fitness cost\n(1-Fitness)") +
  ylab("Trans fitness cost\n(1-Fitness)")
ggsave("Fitnessscatall.pdf", plot = last_plot(), path = figdir, width = 7, height = 7)

nrow(transsizes[transsizes$FitnessJZstrain < transsizes$cisfit,])
nrow(transsizes[transsizes$FitnessJZstrain > transsizes$cisfit,])

length(unique(transsizes[transsizes$trans_med <= transsizes$cisfit,"Focalgene"]))
length(unique(transsizes[transsizes$trans_med >= transsizes$cisfit,"Focalgene"]))

#Overall relationship between fitness and pleiotropy
ggplot(data = cissize, aes(x = log10(cissize), y = FitnessJZstrain)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue2") +
  THEMEMAIN() +
  xlab("Pleiotropy\nNumber of DEGs (Log10)") +
  ylab("Fitness")
ggsave("Fitmod.pdf", plot = last_plot(), path = figdir, width = 7, height = 7)

fitmod <- lm(cissize$FitnessJZstrain ~ log10(cissize$cissize))

#Curious about differences
for (i in 1:nrow(transsizes)) {
  transsizes[i,"Differencemed"] <- transsizes[i,"trans_med"] - transsizes[i,"cisfit"]
  transsizes[i,"Differenceall"] <- transsizes[i,"FitnessJZstrain"] - transsizes[i,"cisfit"]
}

ggplot(data = transsizes, aes(x = Differenceall)) +
  geom_histogram(fill = "grey") +
  geom_vline(xintercept = median(transsizes$Differenceall, na.rm = T), color = "darkred", size = 3) +
  geom_vline(xintercept = 0, size = 3) +
  THEMEMAIN() +
  xlab("trans fitness minus cis fitness") +
  ylab("Count")
ggsave("Fitnesshistall.pdf", plot = last_plot(), path = figdir, width = 7, height = 7)

ggplot(data = transsizes, aes(x = Differencemed)) +
  geom_histogram(fill = "grey") +
  geom_vline(xintercept = median(transsizes$Differencemed, na.rm = T), color = "darkred", size = 3) +
  geom_vline(xintercept = 0, size = 3) +
  THEMEMAIN() +
  xlab("Median trans fitness minus cis fitness") +
  ylab("Count")
ggsave("Fitnesshistmed.pdf", plot = last_plot(), path = figdir, width = 7, height = 7)
