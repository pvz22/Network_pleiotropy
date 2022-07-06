library(pheatmap)
library(ggplot2)
library(dplyr)
library(data.table)
library(gtools)
rm(list = ls())
figdir <- "/Users/petravandezande/Documents/Figures/Projects/Pleiotropy/Friend"
data.dir <- "/Users/petravandezande/Documents/Data/Projects/Pleiodel"
outputdir <- "/Users/petravandezande/Documents/Output/Projects/Pleiotropy/Friend"

THEMEMAIN <- function() {
  theme_bw() +
    theme(axis.text = element_text(size = 25), axis.title = element_text(size = 25), plot.margin = unit(c(2,1,1,1),"cm"), plot.title = element_text(size = 25, hjust = 0.5))
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}
#Skip to line 59 for pleiotropy matrix

mutantsexwtvar <- fread(paste0(data.dir,"/deleteome_all_mutants_ex_wt_var_controls.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)

mutsvswt <- t(mutantsexwtvar[,4:ncol(mutantsexwtvar)])

Mvalues <- mutsvswt[mutsvswt[,1] == "M",]
Mvalues <- Mvalues[,2:ncol(Mvalues)]
Mvalues <- apply(Mvalues, c(1,2), as.numeric)
Mvalues <- as.data.frame(Mvalues, stringsAsFactors = FALSE)
exp.names.vec <- c(rownames(Mvalues))
temp <- strsplit(exp.names.vec, "-del")
temp2 <- lapply(temp, "[[", 1)
mat <- matrix(unlist(temp2), ncol = 1, byrow = TRUE)
mat.df <- as.data.frame(mat, stringsAsFactors = FALSE)
rownames(Mvalues) <- mat.df$V1
colnames(Mvalues) <- mutantsexwtvar$geneSymbol[2:6124]
write.table(Mvalues, paste0(outputdir,"/Mvalues.txt"), sep = "\t")

p.values <- mutsvswt[mutsvswt[,1] == "p_value",]
p.values <- p.values[,2:ncol(p.values)]
p.values <- apply(p.values, c(1,2), as.numeric)
exp.names.vec <- c(rownames(p.values))
temp <- strsplit(exp.names.vec, "-del")
temp2 <- lapply(temp, "[[", 1)
mat <- matrix(unlist(temp2), ncol = 1, byrow = TRUE)
mat.df <- as.data.frame(mat)
rownames(p.values) <- mat.df$V1
p.values <- as.data.frame(p.values, stringsAsFactors = FALSE)
colnames(p.values) <- mutantsexwtvar$geneSymbol[2:6124]
write.table(p.values, paste0(outputdir,"/p.values.txt"), sep = "\t")

Pleiotropy <- matrix(NA, nrow = nrow(Mvalues), ncol = ncol(Mvalues))
Pleiotropy <- ifelse(abs(Mvalues) >= 0.7655347  & pvalues <= 0.05, 1, 0) #Mvalues are actually log2ratios, so FC of 1.7 = 0.7655347
#I end up with 3,000 more than in the Kemmeren publication. This may be due to the eqality or that I didn't filter out 'nonresponsive' experiments.

colnames(Pleiotropy) <- colnames(Mvalues)
rownames(Pleiotropy) <- rownames(Mvalues)
colnames(Pleiotropy) <- tolower(colnames(Pleiotropy))
write.table(Pleiotropy, paste0(outputdir,"/Pleiotropy.txt"), sep = "\t")
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

#Looking at a binary heatmap to see the overall network structure
#Ordering so focal genes will in the same order in rows and columns
formap <- Pleiotropy[rownames(Pleiotropy) %in% colnames(Pleiotropy),colnames(Pleiotropy) %in% rownames(Pleiotropy)]
formap <- formap[sort(rownames(formap)),sort(colnames(formap))]
table(rownames(formap) == colnames(formap))
formap2 <- Pleiotropy[rownames(Pleiotropy) %in% colnames(Pleiotropy),!colnames(Pleiotropy) %in% rownames(Pleiotropy)]
formap2 <- formap2[sort(rownames(formap2)),]
formap3 <- cbind(formap, formap2)
formapmat <- as.matrix(formap3)
plota <- pheatmap(formapmat, border_color = NA, color = c("white","black"), kmeans_k = NA, cellwidth = 1, cellheight = 1, scale = "none", legend = FALSE, show_rownames = 0, show_colnames = 0, cluster_rows = 0, cluster_cols = 0)
ggsave("Amatplot.pdf", plot = plota, path = figdir, height = 18.75, width = 90, limitsize = FALSE)
ggsave("Amatplot.png", plot = plota, path = figdir, height = 18.75, width = 90, limitsize = FALSE)

#Overall distribution of number of trans-regulators for each focal gene
transnum <- data.table(colSums(formap))
max(transnum)
ggplot(data = transnum, aes(x = log10(V1))) +
  geom_histogram(fill = "#336879") +
  THEMEMAIN() +
  xlab("Number of trans-regulators") +
  ylab("Count")
ggsave("transhistlog.pdf", plot = last_plot(), path = figdir, width = 7, height = 7)

ggplot(data = transnum, aes(x = V1)) +
  geom_histogram(fill = "#336879") +
  THEMEMAIN() +
  xlab("Number of trans-regulators") +
  ylab("Count") +
  xlim(0,20)
ggsave("transhistinset.pdf", plot = last_plot(), path = figdir, width = 6, height = 6)

#
Pleiotropy <- Pleiotropy[rownames(Pleiotropy) %in% colnames(Pleiotropy),] #We can only look at focal genes that we can ID trans regulators for
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
cissize <- data.frame(cissize)
write.table(transsizes, paste0(outputdir,"/transsizeslim.txt"), sep = "\t")
write.table(cissize, paste0(outputdir,"/cissizelim.txt"), sep = "\t")
transsizes$Focalgene <- as.character(transsizes$Focalgene)
transsizes$Strain <- as.character(transsizes$Strain)
#Removing self-loops as trans-regulatory perturbations - 
#also checking that all have self-loops to begin with!!
table(unique(transsizes$Focalgene) %in% unique(transsizes$Strain))
transsizes <- transsizes[which(transsizes$Strain != transsizes$Focalgene),]

#Total distribution of cis-acting pleiotropy
min(cissize[!is.na(cissize)])
ggplot(data = cissize, aes(x = cissize)) +
  geom_histogram(fill = "seagreen4") +
  THEMEMAIN() +
  xlab("Cis pleiotropy\n(Number of DEGs)") +
  ylab("Count")
ggsave("Histcis.pdf", plot = last_plot(), path = figdir, width = 7, height = 7)

ggplot(data = cissize, aes(x = cissize)) +
  geom_histogram(fill = "seagreen4") +
  THEMEMAIN() +
  xlab("Cis pleiotropy\n(Number of DEGs)") +
  ylab("Count") +
  xlim(-5,100)
ggsave("Histcisinset.pdf", plot = last_plot(), path = figdir, width = 6, height = 6)


#Looking at an example of Ptr2
ggplot(data = transsizes[transsizes$Focalgene == "ptr2",], aes(x = Size)) +
  geom_histogram(fill = "darkorchid4", alpha = 0.8) +
  geom_vline(xintercept = median(transsizes[transsizes$Focalgene == "ptr2","Size"]), color = "black", size = 3) +
  geom_vline(xintercept = cissize["ptr2",], color = "seagreen4", size = 3) +
  THEMEMAIN() +
  xlab("Trans pleiotropy\n(Number of DEGs)") +
  ylab("Count") +
  ggtitle("PTR2")
ggsave("ptr2.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)

#Creating a single value for all trans-regs for each focal gene
trans_avg <- transsizes %>%
  group_by(Focalgene) %>%
  summarize(avg = mean(Size))

trans_med <- transsizes %>%
  group_by(Focalgene) %>%
  summarize(med = median(Size))


#Scatterplot of cis vs trans
for (i in 1:nrow(transsizes)) {
  transsizes[i,"cissize"] <- cissize[rownames(cissize) == transsizes[i,"Focalgene"],"cissize"]
  transsizes[i,"trans_med"] <- trans_med[trans_med$Focalgene == transsizes[i,"Focalgene"],"med"]
  transsizes[i,"trans_avg"] <- trans_avg[trans_avg$Focalgene == transsizes[i,"Focalgene"],"avg"]
}
#Reading this out to save it
write.table(transsizes, paste0(outputdir,"/transsizeslim.txt"), sep = "\t")

ggplot(data = transsizes, aes(x = log10(cissize), y = log10(trans_med))) +
  geom_point(color = "black") +
  geom_abline(slope = 1, intercept = 0, color = "seagreen4", size = 2) +
  THEMEMAIN() +
  xlab("Cis pleiotropy\n(Log10)") +
  ylab("Median trans pleiotropy\n(Log10)") +
  xlim(0,3.1) +
  ylim(0,3.1)
ggsave("Scatterplotcisvstransmedlog.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)

#All trans-regulators separately rather than the median
ggplot(data = transsizes, aes(x = log10(cissize), y = log10(Size))) +
  geom_point(color = "darkorchid4") +
  geom_abline(slope = 1, intercept = 0, color = "seagreen4", size = 2) +
  THEMEMAIN() +
  xlab("Cis pleiotropy\n(Log10)") +
  ylab("Trans pleiotropy\n(Log10)")
ggsave("Scatterplotcisvstranslog.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)

#Looking at a histogram of all the differences, or essentially residuals from the previous plots
trans_avg <- as.data.frame(trans_avg)
for (i in 1:nrow(trans_avg)) {
  trans_avg[i,"Difference"] <- trans_avg[i,"avg"] - cissize[as.character(trans_avg[i,"Focalgene"]),"cissize"]
}
trans_med <- as.data.frame(trans_med)
for (i in 1:nrow(trans_med)) {
  trans_med[i,"Difference"] <- trans_med[i,"med"] - cissize[as.character(trans_med[i,"Focalgene"]),"cissize"]
}

# ggplot(data = trans_avg, aes(x = Difference)) +
#   geom_histogram(aes(x = Difference), fill = "grey", alpha = 0.7) +
#   THEMEMAIN() +
#   xlab("Difference in number \nof DE genes in trans and cis") +
#   geom_vline(xintercept = mean(trans_avg$Difference), color = "#336879", size = 2) +
#   geom_vline(xintercept = 0, color = "black", size = 2) +
#   ylab("Density")
# ggsave("Densitytransminuscis.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)

ggplot(data = trans_med, aes(x = Difference)) +
  geom_histogram(aes(x = Difference), fill = "grey") +
  THEMEMAIN() +
  xlab("Median trans minus cis") +
  ylab("Count") +
  geom_vline(xintercept = mean(trans_med$Difference), color = "darkred", size = 2) +
  geom_vline(xintercept = 0, color = "black", size = 2)
ggsave("Histtransminuscismed.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)
t.test(trans_med$Difference, alternative = "greater")$p.value #p-value: < 2.2e-16
nrow(trans_med[trans_med$Difference < 0,]) #26 out of 804, 3%

#A histogram of the differences for each individual trans-regulator and cis rather than sum stats
transsizes$Difference <- transsizes$Size - transsizes$cissize
ggplot(data = transsizes, aes(x = Difference)) +
  geom_histogram(fill = "darkorchid4", alpha = 0.7) +
  THEMEMAIN() +
  xlab("Each trans minus cis") +
  ylab("Count") +
  geom_vline(xintercept = mean(transsizes$Difference), color = "darkred", size = 2) +
  geom_vline(xintercept = 0, color = "black", size = 2)
ggsave("Histtransminuscis.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)
t.test(transsizes$Difference, alternative = "greater")$p.value 
nrow(transsizes[transsizes$Difference < 0,]) #365 out of 8392, 4%

#Power-law relationship for the out-degree distribution
#calculating p(K) for each K value:
Kvalues <- rowSums(Pleiotropy)
Ktable <- table(Kvalues)
Ktable <- data.frame(Ktable)
Ktable$Kvalues <- as.numeric(as.character(Ktable$Kvalues))
Ktable$pK <- Ktable$Freq/sum(Ktable$Freq)
#Getting rid of the four zeros that mess with the log
Ktable <- Ktable[2:nrow(Ktable),]
ggplot(data = Ktable, aes(x = log10(Kvalues), y = log10(pK))) +
  geom_point(size = 3) +
  THEMEMAIN() +
  xlab("Out-Degree\n(Log10)") +
  ylab("Probably of Out-Degree\n(Log10)") +
  geom_smooth(method = "lm")
  #geom_abline(slope = -0.75, intercept = -1.3)
ggsave("kpkplot.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)

fit <- lm(log10(Ktable$pK) ~ log10(Ktable$Kvalues))
Ktable$approx <- Ktable$Kvalues^-0.75
plot(log10(Ktable$Kvalues), log10(Ktable$pK))
abline(fit)
points(log10(Ktable$Kvalues), log10(Ktable$approx) - 1, col = "red")



#Histogram of trans sizes all together
ggplot(data = transsizes, aes(x = Size)) +
  geom_histogram(aes(x = Size), fill = "darkorchid4") +
  THEMEMAIN() +
  xlab("Trans pleiotropy\n(Number of DEGs)") +
  ylab("Count")
ggsave("Histtrans.pdf", plot = last_plot(), path = figdir, width = 7, height = 7)



#Looking at effect sizes for all focal genes in cis and in trans
for (i in 1:nrow(transsizes)) {
  transsizes[i,"Effectsize"] <- Mvalues[rownames(Mvalues) == transsizes[i,"Strain"],tolower(colnames(Mvalues)) == transsizes[i,"Focalgene"]]
}
for (i in 1:nrow(transsizes)) {
  transsizes[i,"Effectsizecis"] <- Mvalues[rownames(Mvalues) == transsizes[i,"Focalgene"],tolower(colnames(Mvalues)) == transsizes[i,"Focalgene"]]
}
ggplot(data = transsizeslim[transsizeslim$Effectsize < transsizeslim$Effectsizecis,], aes(x = Focalgene, y = Effectsize)) +
  geom_point(color = "#336879", alpha = 0.7) +
  geom_point(aes(x = Focalgene, y = Effectsizecis), color = "#f47a60", alpha = 0.7) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_blank(), legend.position = "none",axis.title = element_text(size = 25), plot.margin = unit(c(2,1,1,1),"cm")) +
  ylab("Focal Gene Expression\n(Log2 Fold Change)") +
  xlab("Focal Gene")
ggsave("Cisvstranseffectweird.pdf", plot = last_plot(), path = figdir, width = 8, height = 7)

# ggplot(data = transsizes[1:1000,], aes(x = Focalgene, y = Effectsize)) +
#   geom_point(color = "#336879") +
#   geom_point(aes(x = Focalgene, y = Effectsizecis), color = "#f47a60")
# ggsave("Cisvstranseffect.pdf", plot = last_plot(), path = figdir, width = 4, height = 3)

ggplot(data = transsizeslim, aes(x = Focalgene, y = Effectsize)) +
  geom_point(color = "#336879", alpha = 0.7) +
  geom_point(aes(x = Focalgene, y = Effectsizecis), color = "#f47a60", alpha = 0.7) +
  theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_blank(), legend.position = "none",axis.title = element_text(size = 25), plot.margin = unit(c(2,1,1,1),"cm")) +
  ylab("Focal Gene Expression\n(Log2 Fold Change)") +
  xlab("Focal Gene")
ggsave("Cisvstranseffectall.pdf", plot = last_plot(), path = figdir, width = 8, height = 7)


length(unique(transsizes[transsizes$Effectsize < transsizes$Effectsizecis,"Focalgene"]))
length(unique(transsizes[,"Focalgene"]))

length(unique(transsizes[transsizes$Effectsizecis <= -0.76,"Focalgene"]))

transsizeslim <- transsizes[transsizes$Effectsizecis <= -0.76,]

table(unique(transsizes[abs(transsizes$Effectsizecis) <= 0.76,"Focalgene"]) %in% unique(transsizes[transsizes$Effectsize < transsizes$Effectsizecis,"Focalgene"]))

