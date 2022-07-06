#Permutations scaled up
library(dplyr)
library(cowplot)
library(ggplot2)
library(gtools)
outputdir <- "/Users/petravandezande/Documents/Output/Projects/Pleiotropy/Friend"
figdir <- "/Users/petravandezande/Documents/Figures/Projects/Pleiotropy/Friend"

THEMEMAIN <- function() {
  theme_bw() +
    theme(axis.text = element_text(size = 25), axis.title = element_text(size = 25), plot.margin = unit(c(2,1,1,1),"cm"), plot.title = element_text(size = 25, hjust = 0.5))
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

Pleiotropy <- read.table(paste0(outputdir,"/Pleiotropy.txt"), sep = "\t")

#Removing focal genes that do not show a significant decrease in expression upon deletion
noreducvec <- c()
for (i in 1:nrow(Pleiotropy)) {
  if(rownames(Pleiotropy)[i] %in% colnames(Pleiotropy)) {
    if(Pleiotropy[i,colnames(Pleiotropy) == rownames(Pleiotropy)[i]] == 0) {
      noreducvec <- c(noreducvec,i)
    } 
  }
}
#Removes 127 genes
Pleiotropy <- Pleiotropy[-noreducvec,]
Pleiotropyori <- Pleiotropy[rownames(Pleiotropy) %in% colnames(Pleiotropy),] #We can only look at focal genes that we can ID trans regulators for

#Now generating permutations that maintain degree distribution (permuting rows)
perm1medians <- c()
b <- 1
while (b < 101) {
  print(b)
  Pleiotropyt <- t(Pleiotropyori)
  for (i in 1:ncol(Pleiotropyt)) {
    Pleiotropyt[,i] <- gtools::permute(Pleiotropyt[,i])
  }
  Pleiotropy <- t(Pleiotropyt)
  cissize <- c()
  transsizes <- data.frame()
  allsizes <- data.frame("Strain" = rownames(Pleiotropy), "Size" = rowSums(Pleiotropy) - 2) #Trans size includes focal gene, and trans regulator itself, so subtract 2
  cissize <- rowSums(Pleiotropy) - 1 #cissize includes the focal gene, so subtract 1
  for (i in 1:nrow(Pleiotropy)) {
    focalgene <- rownames(Pleiotropy)[i]
    if(sum(Pleiotropy[,focalgene]) > 1) { #To eliminate where it is only the self loop as the trans regulator
      transregs <- rownames(Pleiotropy[Pleiotropy[,focalgene] == 1,])
      trans <- allsizes[rownames(allsizes) %in% transregs,]
      trans$Focalgene <- rep(focalgene, nrow(trans))
      transsizes <- rbind(transsizes, trans)
    }
  }
  cissize <- data.frame("Cissize" = cissize)
  transsizes$Focalgene <- as.character(transsizes$Focalgene)
  transsizes$Strain <- as.character(transsizes$Strain)
  #Removing self-loops as trans-regulatory perturbations
  transsizes <- transsizes[which(transsizes$Strain != transsizes$Focalgene),]
  trans_med <- transsizes %>%
    group_by(Focalgene) %>%
    summarize(med = median(Size))
  trans_med <- as.data.frame(trans_med)
  for (i in 1:nrow(trans_med)) {
    trans_med[i,"Difference"] <- trans_med[i,"med"] - cissize[as.character(trans_med[i,"Focalgene"]),"Cissize"]
  }
  perm1medians[b] <- median(trans_med$Difference)
  b <- b + 1
}

ggplot(data = as.data.frame(perm1medians), aes(x = perm1medians)) +
  geom_histogram(fill = "darkgrey", bins = 100) +
  THEMEMAIN() +
  geom_vline(xintercept = 0, color = "black", size = 3) +
  ylab("Count") +
  xlab("Median difference\n in cis and trans pleiotropy\n(100 permutations)")
ggsave("perm1meds.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)

#Permutaing all edges, not maintaining degree distribution
perm2medians <- c()
b <- 1
while (b < 101) {
  print(b)
  Pleiotropyt <- t(Pleiotropyori)
  for (i in 1:ncol(Pleiotropyt)) {
    Pleiotropyt[,i] <- gtools::permute(Pleiotropyt[,i])
  }
  Pleiotropy <- t(Pleiotropyt)
  for (i in 1:ncol(Pleiotropy)) {
    Pleiotropy[,i] <- gtools::permute(Pleiotropy[,i])
  }
  cissize <- c()
  transsizes <- data.frame()
  allsizes <- data.frame("Strain" = rownames(Pleiotropy), "Size" = rowSums(Pleiotropy) - 2) #Trans size includes focal gene, and trans regulator itself, so subtract 2
  cissize <- rowSums(Pleiotropy) - 1 #cissize includes the focal gene, so subtract 1
  for (i in 1:nrow(Pleiotropy)) {
    focalgene <- rownames(Pleiotropy)[i]
    if(sum(Pleiotropy[,focalgene]) > 1) { #To eliminate where it is only the self loop as the trans regulator
      transregs <- rownames(Pleiotropy[Pleiotropy[,focalgene] == 1,])
      trans <- allsizes[rownames(allsizes) %in% transregs,]
      trans$Focalgene <- rep(focalgene, nrow(trans))
      transsizes <- rbind(transsizes, trans)
    }
  }
  cissize <- data.frame("Cissize" = cissize)
  transsizes$Focalgene <- as.character(transsizes$Focalgene)
  transsizes$Strain <- as.character(transsizes$Strain)
  #Removing self-loops as trans-regulatory perturbations
  transsizes <- transsizes[which(transsizes$Strain != transsizes$Focalgene),]
  trans_med <- transsizes %>%
    group_by(Focalgene) %>%
    summarize(med = median(Size))
  trans_med <- as.data.frame(trans_med)
  for (i in 1:nrow(trans_med)) {
    trans_med[i,"Difference"] <- trans_med[i,"med"] - cissize[as.character(trans_med[i,"Focalgene"]),"Cissize"]
  }
  perm2medians[b] <- median(trans_med$Difference)
  b <- b + 1
}

#Plots for one of these permutations
ggplot(data = trans_med, aes(x = Difference)) +
  geom_histogram(fill = "darkgrey") +
  THEMEMAIN() +
  geom_vline(xintercept = 0, color = "black", size = 2) +
  geom_vline(xintercept = median(trans_med$Difference), color = "darkred", size = 2) +
  ylab("Count") +
  xlab("Trans pleiotropy minus Cis")
ggsave("Transminuscisperm2.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)

ggplot(data = cissize, aes(x = Cissize)) +
  geom_histogram(fill = "seagreen4") +
  THEMEMAIN() +
  xlab("Cis pleiotropy\n(Number of DEGs)") +
  ylab("Count")
ggsave("Histcisperm2.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)

ggplot(data = transsizes, aes(x = Size)) +
  geom_histogram(fill = "darkorchid4") +
  THEMEMAIN() +
  xlab("Trans pleiotropy\n(Number of DEGs)") +
  ylab("Count")
ggsave("Histtransperm2.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)

#Plot of medians for all 100 of these permutations
ggplot(data = as.data.frame(perm2medians), aes(x = perm2medians)) +
  geom_histogram(fill = "darkgrey") +
  THEMEMAIN() +
  ylab("Count") +
  xlab("Median difference\n in cis and trans\n(100 permutations)")
ggsave("perm2meds.pdf", plot = last_plot(), path = figdir, width = 4, height = 4)

