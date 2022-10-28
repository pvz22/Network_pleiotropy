#Scripts used to analyze data from Kemmeren et al (2014) and generate figures for Vande Zande et al 2022, MBE
#Libraries used
library(pheatmap)
library(ggplot2)
library(dplyr)
library(data.table)
library(gtools)
#Setting up environmental variables
rm(list = ls())
figdir <- "/Users/petravandezande/Documents/Figures/Projects/Pleiotropy/Friend/Revision"
data.dir <- "/Users/petravandezande/Documents/Data/Projects/Pleiodel"
outputdir <- "/Users/petravandezande/Documents/Output/Projects/Pleiotropy/Friend"
#GGplot theme
THEMEMAIN <- function() {
  theme_bw() +
    theme(axis.text = element_text(size = 25), axis.title = element_text(size = 25), plot.margin = unit(c(2,1,1,1),"cm"), plot.title = element_text(size = 25, hjust = 0.5))
}
#Reading in the data from Kemmeren et al (2014) and extracting M-values and p-values
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
Mvalues <- read.table(paste0(outputdir,"/Mvalues.txt"), sep = "\t", header = 1)
colnames(Mvalues) <- tolower(colnames(Mvalues))

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
p.values <- read.table(paste0(outputdir,"/p.values.txt"), sep = "\t", header = 1)
colnames(p.values) <- tolower(colnames(p.values))

#Generating a binary matrix of DEGs in each deletion strain ()
#First checking that all things are positioned correctly
table(rownames(Mvalues) == rownames(p.values))
table(colnames(Mvalues) == colnames(p.values))

Pleiotropy <- matrix(NA, nrow = nrow(Mvalues), ncol = ncol(Mvalues))
Pleiotropy <- ifelse(abs(Mvalues) >= 0.7655347  & p.values <= 0.05, 1, 0) #Mvalues are log2ratios, so FC of 1.7 = 0.7655347

colnames(Pleiotropy) <- colnames(Mvalues)
rownames(Pleiotropy) <- rownames(Mvalues)
write.table(Pleiotropy, paste0(outputdir,"/Pleiotropy.txt"), sep = "\t")
Pleiotropy <- read.table(paste0(outputdir,"/Pleiotropy.txt"), sep = "\t")

#Removing focal genes that do not show a significant decrease in their own expression upon deletion
noreducvec <- c()
for (i in 1:nrow(Pleiotropy)) {
  if(rownames(Pleiotropy)[i] %in% colnames(Pleiotropy)) {
    if(Pleiotropy[i,colnames(Pleiotropy) == rownames(Pleiotropy)[i]] == 0) {
      noreducvec <- c(noreducvec,i)
    }
  }
}
for (i in 1:nrow(Mvalues)) {
  if(rownames(Mvalues)[i] %in% colnames(Mvalues)) {
    if(Mvalues[i,colnames(Mvalues) == rownames(Mvalues)[i]] > 0) {
      noreducvec <- c(noreducvec,i)
    }
  }
}
#Removes 129 genes
Pleiotropy <- Pleiotropy[-noreducvec,]

#Removing deletions whose expression were not assayed by the microarray (removes 82 genes, 3 non-deletion measures)
Pleiotropy <- Pleiotropy[rownames(Pleiotropy) %in% colnames(Pleiotropy),]
write.table(Pleiotropy,paste0(outputdir,"/TableS2.txt"), sep = "\t", row.names = TRUE)

#Calculating number of DEGs for cis and trans-regulators for each focal gene
transsizes <- data.frame()
allsizes <- data.frame("Strain" = rownames(Pleiotropy), "Size" = rowSums(Pleiotropy) - 2) #Trans size includes focal gene, and trans regulator itself, so subtract 2
cissize <- rowSums(Pleiotropy) - 1 #cissize includes the focal gene, so subtract 1
translist <- list()
for (i in 1:nrow(Pleiotropy)) {
  focalgene <- rownames(Pleiotropy)[i]
  transregsvec <- rownames(Pleiotropy[Pleiotropy[,focalgene] == 1,])
  trans <- allsizes[rownames(allsizes) %in% transregsvec,]
  trans$Focalgene <- rep(focalgene, nrow(trans))
  transsizes <- rbind(transsizes, trans)
  transregs <- list(rownames(Pleiotropy[Pleiotropy[,focalgene] == 1,]))
  names(transregs) <- focalgene
  translist <- append(translist, transregs)
}
transsizes$Focalgene <- as.character(transsizes$Focalgene)
transsizes$Strain <- as.character(transsizes$Strain)

#double checking that all focal genes have a 'self-loop' or are DE in their own deletion
table(unique(transsizes$Focalgene) %in% unique(transsizes$Strain))
#Removing instances where the gene is only DE in its own deletion (aka, has no other trans-regulators)
transsizes <- transsizes[which(transsizes$Strain != transsizes$Focalgene),]
length(unique(transsizes$Focalgene)) #748, removes 525

#Reading out the pleiotropy adjacency matrix containing only the 748 focal genes
Pleiotropyfocalgenes <- Pleiotropy[rownames(Pleiotropy) %in% unique(transsizes$Focalgene),]
write.table(Pleiotropyfocalgenes, paste0(outputdir,"/Pleiotropyfocalgenes.txt"), sep = "\t")

#Total distribution of cis-acting pleiotropy
cissize <- cissize[unique(transsizes$Focalgene)] #Only looking at the distribution for focal genes that we have ALL the information for
cissize <- data.frame(cissize = cissize)
write.table(cissize, paste0(outputdir,"/cissizelim.txt"), sep = "\t")
min(cissize[!is.na(cissize$cissize),"cissize"])
max(cissize[!is.na(cissize$cissize),"cissize"])

ggplot(data = cissize, aes(x = cissize)) +
  geom_histogram(fill = "seagreen4") +
  THEMEMAIN() +
  xlab("Cis-regulatory pleiotropy\n(Number of DEGs)") +
  ylab("Count")
ggsave("Histcis.pdf", plot = last_plot(), path = figdir, width = 7, height = 7)

ggplot(data = cissize, aes(x = cissize)) +
  geom_histogram(fill = "seagreen4") +
  THEMEMAIN() +
  xlab("Cis-regulatory pleiotropy\n(Number of DEGs)") +
  ylab("Count") +
  xlim(-5,100)
ggsave("Histcisinset.pdf", plot = last_plot(), path = figdir, width = 6, height = 6)

#Overall distribution of number of trans-regulators for each focal gene
transnumvec <-c(table(transsizes$Focalgene))
max(transnumvec)
min(transnumvec)
transnum <- data.frame(Transnum = transnumvec)

ggplot(data = transnum, aes(x = Transnum)) +
  geom_histogram(fill = "#336879", binwidth = 10) +
  THEMEMAIN() +
  xlab("Number of trans-regulators") +
  ylab("Count")
ggsave("transhist.pdf", plot = last_plot(), path = figdir, width = 7, height = 7)

ggplot(data = transnum, aes(x = Transnum)) +
  geom_histogram(fill = "#336879", binwidth = 1) +
  THEMEMAIN() +
  xlab("Number of trans-regulators") +
  ylab("Count") +
  xlim(0,20)
ggsave("transhistinset.pdf", plot = last_plot(), path = figdir, width = 6, height = 6)

#Looking at the relationship between indegree and outdegree
cissize$Gene <- rownames(cissize)
transnum$Gene <- rownames(transnum)
Degrees <- left_join(transnum, cissize, by = "Gene")

ggplot(data = Degrees, aes(x = cissize, y = Transnum)) +
  geom_point(size = 5) +
  THEMEMAIN() +
  xlab("Cis-regulatory pleiotropy\n (Number of DEGs)") +
  ylab("Number of trans-regulators")
ggsave("cispleiovstransnum.pdf", plot = last_plot(), path = figdir, width = 7, height = 7)


#Looking at an example of Ptr2
ggplot(data = transsizes[transsizes$Focalgene == "ptr2",], aes(x = Size)) +
  geom_histogram(fill = "darkorchid4", alpha = 0.8) +
  geom_vline(xintercept = median(transsizes[transsizes$Focalgene == "ptr2","Size"]), color = "black", size = 3) +
  geom_vline(xintercept = cissize["ptr2","cissize"], color = "seagreen4", size = 3) +
  THEMEMAIN() +
  xlab("Trans-regulatory pleiotropy\n(Number of DEGs)") +
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

ggplot(data = transsizes, aes(x = log10(cissize), y = log10(trans_med))) +
  geom_point(color = "black") +
  geom_abline(slope = 1, intercept = 0, color = "seagreen4", size = 2) +
  THEMEMAIN() +
  xlab("Cis-regulatory pleiotropy\n(Log10 Number of DEGs)") +
  ylab("Median trans-regulatory pleiotropy\n(Log10 Number of DEGs)") +
  xlim(0,3.1) +
  ylim(0,3.1)
ggsave("Scatterplotcisvstransmedlog.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)


#Looking at a histogram of all the differences, or essentially residuals from the previous plots
trans_avg <- as.data.frame(trans_avg)
for (i in 1:nrow(trans_avg)) {
  trans_avg[i,"Difference"] <- trans_avg[i,"avg"] - cissize[as.character(trans_avg[i,"Focalgene"]),"cissize"]
}
trans_med <- as.data.frame(trans_med)
for (i in 1:nrow(trans_med)) {
  trans_med[i,"Difference"] <- trans_med[i,"med"] - cissize[as.character(trans_med[i,"Focalgene"]),"cissize"]
}

ggplot(data = trans_med, aes(x = Difference)) +
  geom_histogram(aes(x = Difference), fill = "grey") +
  THEMEMAIN() +
  xlab("Median trans-regulatory pleiotropy\n minus cis-regulatory pleiotropy\n (Log10 Number of DEGs)") +
  ylab("Count") +
  geom_vline(xintercept = mean(trans_med$Difference), color = "darkred", size = 2) +
  geom_vline(xintercept = 0, color = "black", size = 2)
ggsave("Histtransminuscismed.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)
t.test(trans_med$Difference, alternative = "greater")$p.value #p-value: 3.0e-165
nrow(trans_med[trans_med$Difference > 0,]) #718 out of 748, ~96%


#All trans-regulators separately rather than the median
ggplot(data = transsizes, aes(x = log10(cissize), y = log10(Size))) +
  geom_point(color = "darkorchid4") +
  geom_abline(slope = 1, intercept = 0, color = "seagreen4", size = 2) +
  THEMEMAIN() +
  xlab("Cis-regulatory pleiotropy\n(Log10 Number of DEGs)") +
  ylab("Trans-regulatory pleiotropy\n(Log10 Number of DEGs)")
ggsave("Scatterplotcisvstranslog.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)
nrow(transsizes[transsizes$Size > transsizes$cissize,]) #6785 of 7162, 95%

#A histogram of the differences for each individual trans-regulator and cis rather than sum stats
transsizes$Difference <- transsizes$Size - transsizes$cissize
ggplot(data = transsizes, aes(x = Difference)) +
  geom_histogram(fill = "darkorchid4", alpha = 0.7) +
  THEMEMAIN() +
  xlab("Each trans-regulatory pleiotropy\n minus cis-regulatory pleiotropy\n (Number of DEGs)") +
  ylab("Count") +
  geom_vline(xintercept = mean(transsizes$Difference), color = "darkred", size = 2) +
  geom_vline(xintercept = 0, color = "black", size = 2)
ggsave("Histtransminuscis.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)
t.test(transsizes$Difference, alternative = "greater")$p.value #0
nrow(transsizes[transsizes$Difference <= 0,]) #377 of 7162

#Looking at all of this with euclidean distances instead of number of DEGs
#First, just making a dataframe with all the euclidean distances for all deletions.
#I will calculate the euclidean distances first without eliminating any genes, straight from the Mvalues matrix
Mvaluessq <- Mvalues^2
Eudists <- rowSums(Mvaluessq)
Eudists <- sqrt(Eudists)
#Now I just have to add these values onto the transsizes dataframe
for (i in 1:nrow(transsizes)) {
  transsizes[i,"Eudisttrans"] <- Eudists[transsizes[i,"Strain"]]
  transsizes[i,"Eudistcis"] <- Eudists[transsizes[i,"Focalgene"]]
}
#Comparing all pairwise euclidean distances
ggplot(data = transsizes, aes(x = Eudistcis, y = Eudisttrans)) +
  geom_point(color = "darkorchid4") +
  geom_abline(slope = 1, intercept = 0, color = "seagreen4", size = 2) +
  THEMEMAIN() +
  xlab("Cis-regulatory pleiotropy\n(Euclidean distance)") +
  ylab("Trans-regulatory pleiotropy\n(Euclidean distance)")
ggsave("Scatterplotcisvstranseudist.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)
nrow(transsizes[transsizes$Eudisttrans <= transsizes$Eudistcis,])#377 of 7162, 5%
#NOTE - the summary stat (median and average) comparisons are below

#Just the distribution of cis pleiotropy as euclidean distances
for (i in 1:nrow(cissize)) {
  cissize[i,"Eudist"] <- Eudists[cissize[i,"Gene"]]
}

ggplot(data = cissize, aes(x = Eudist)) +
  geom_histogram(fill = "seagreen4") +
  THEMEMAIN() +
  xlab("Cis-regulatory pleiotropy\n(Euclidean distance)") +
  ylab("Count")
ggsave("Histciseudist.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)

#Now doing the summaries of trans values
trans_avg_eud <- transsizes %>%
  group_by(Focalgene) %>%
  summarize(avgeud = mean(Eudisttrans))

trans_med_eud <- transsizes %>%
  group_by(Focalgene) %>%
  summarize(medeud = median(Eudisttrans))


#Putting back into the original table
for (i in 1:nrow(transsizes)) {
  transsizes[i,"trans_eud_med"] <- trans_med_eud[trans_med_eud$Focalgene == transsizes[i,"Focalgene"],"medeud"]
  transsizes[i,"trans_eud_avg"] <- trans_avg_eud[trans_avg_eud$Focalgene == transsizes[i,"Focalgene"],"avgeud"]
}

#Plotting median euclidean distances for cis and trans
ggplot(data = transsizes, aes(x = Eudistcis, y = trans_eud_med)) +
  geom_point(color = "black") +
  geom_abline(slope = 1, intercept = 0, color = "seagreen4", size = 2) +
  THEMEMAIN() +
  xlab("Cis-regulatory pleiotropy\n(Euclidean distance)") +
  ylab("Median trans-regulatory pleiotropy\n(Euclidean distance)")
ggsave("Scatterplotcisvstranseudistmed.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)

trans_med_eud <- as.data.frame(trans_med_eud)
for (i in 1:nrow(trans_med_eud)) {
  trans_med_eud[i,"Difference"] <- trans_med_eud[i,"medeud"] - cissize[as.character(trans_med_eud[i,"Focalgene"]),"Eudist"]
}
nrow(trans_med_eud[trans_med_eud$Difference <= 0,]) #32 out of 748, 4%

#Is trans still bigger if we remove cis-effects?
#Identifying the downstream effects of each cis mutant that are also present in trans mutants
downstreamlist <- list()
for (i in 1:nrow(Pleiotropy)) {
  downstreamlist[[i]] <- c(colnames(Pleiotropy[,Pleiotropy[i,] == 1]))
}
names(downstreamlist) <- rownames(Pleiotropy)
for (i in 1:nrow(transsizes)) {
  transsizes[i,"Numnested"] <- table(downstreamlist[[transsizes[i,"Strain"]]] %in% downstreamlist[[transsizes[i,"Focalgene"]]])["TRUE"]
}
#These numbers of nested effects include the focal gene (so we have to subtract 1). Also, if there was only the focal gene nested,
#This code will input an "NA", which is fine because there are no nested effects to subtract in that case
transsizes$Numnested <- ifelse(is.na(transsizes$Numnested),1,transsizes$Numnested)
trans_nest_med <- transsizes %>%
  mutate(Sizeminusnest = Size - (Numnested - 1)) %>%
  group_by(Focalgene) %>%
  summarise(mednestsub = median(Sizeminusnest))

transsizes <- data.frame(left_join(transsizes, trans_nest_med, by = "Focalgene"))

ggplot(data = transsizes, aes(x = log10(cissize), y = log10(mednestsub))) +
  geom_point(color = "black") +
  geom_abline(slope = 1, intercept = 0, color = "seagreen4", size = 2) +
  THEMEMAIN() +
  xlab("Cis-regulatory pleiotropy\nLog10(Number of DEGs)") +
  ylab("Median trans-regulatory pleiotropy\nin parallel\nLog10(Number of DEGs - Nested effects)")
ggsave("Scatterplotcisvstransmednonest.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)
length(unique(transsizes[transsizes$cissize <= transsizes$mednestsub,"Focalgene"])) #712 out of 748
nrow(transsizes[transsizes$cissize >= (transsizes$Size - (transsizes$Numnested - 1)),]) #413 out of 7162


#Looking at plots for just the trans-regulatory factors that reduce focal gene expression
for (i in 1:nrow(transsizes)) {
  transsizes[i,"Effectsize"] <- Mvalues[transsizes[i,"Strain"],transsizes[i,"Focalgene"]]
}

downreg <- transsizes[transsizes$Effectsize < 0,]

ggplot(data = downreg, aes(x = log10(cissize), y = log10(Size - (Numnested - 1)))) +
  geom_point(color = "darkorchid4") +
  geom_abline(slope = 1, intercept = 0, color = "seagreen4", size = 2) +
  THEMEMAIN() +
  xlab("Cis-regulatory pleiotropy\nLog10(Number of DEGs)") +
  ylab("Down-regulator trans-regulatory pleiotropy\n(Number of DEGs - nested effects)")
ggsave("Scatterplotcisvstransnestdown.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)
nrow(downreg[downreg$cissize >= (downreg$Size - (downreg$Numnested - 1)),]) #192 out of 2257

ggplot(data = downreg, aes(x = Eudistcis, y = Eudisttrans)) +
  geom_point(color = "darkorchid4") +
  geom_abline(slope = 1, intercept = 0, color = "seagreen4", size = 2) +
  THEMEMAIN() +
  xlab("Cis-regulatory pleiotropy\n(Euclidean distance)") +
  ylab("Down-regulator trans-regulatory pleiotropy\n(Euclidean distance)")
ggsave("Scatterplotcisvstranseuddown.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)

ggplot(data = downreg, aes(x = log10(cissize), y = log10(Size))) +
  geom_point(color = "darkorchid4") +
  geom_abline(slope = 1, intercept = 0, color = "seagreen4", size = 2) +
  THEMEMAIN() +
  xlab("Cis-regulatory pleiotropy\n(Log10 Number of DEGs)") +
  ylab("Down-regulator trans-regulatory pleiotropy\n(Log10 Number of DEGs)")
ggsave("Scatterplotcisvstransdown.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)

#Now looking at euclidean distances with the slow growth signature removed
sigsub <- fread(paste0(data.dir,"/deleteome_all_mutants_svd_transformed.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
genenames <- sigsub[,2]
sigsub <- sigsub[,3:ncol(sigsub)]
sigsub <- apply(sigsub, c(1,2), as.numeric)
sigsub <- as.data.frame(t(sigsub), stringsAsFactors = FALSE)
exp.names.vec <- c(rownames(sigsub))
temp <- strsplit(exp.names.vec, ".del")
temp2 <- lapply(temp, "[[", 1)
mat <- matrix(unlist(temp2), ncol = 1, byrow = TRUE)
mat.df <- as.data.frame(mat, stringsAsFactors = FALSE)
rownames(sigsub) <- mat.df$V1
colnames(sigsub) <- genenames$commonName
colnames(sigsub) <- tolower(colnames(sigsub))

sigsubsq <- sigsub^2
Eudistssub <- rowSums(sigsubsq)
Eudistssub <- sqrt(Eudistssub)
#Now I just have to add these values onto the transsizes dataframe
for (i in 1:nrow(transsizes)) {
  transsizes[i,"Eudisttranssub"] <- Eudistssub[transsizes[i,"Strain"]]
  transsizes[i,"Eudistcissub"] <- Eudistssub[transsizes[i,"Focalgene"]]
}

#Okay, now comparing these - looking at all pairwise first
ggplot(data = transsizes, aes(x = Eudistcissub, y = Eudisttranssub)) +
  geom_point(color = "darkorchid4") +
  geom_abline(slope = 1, intercept = 0, color = "seagreen4", size = 2) +
  THEMEMAIN() +
  xlab("Cis-regulatory pleiotropy\n(Euclidean distance)\nSlow growth removed") +
  ylab("Trans-regulatory pleiotropy\n(Euclidean distance)\nSlow growth removed")
ggsave("Scatterplotcisvstranseudistsub.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)

#Now doing the summaries of trans values
trans_avg2 <- transsizes %>%
  group_by(Focalgene) %>%
  summarize(avgeudsub = mean(Eudisttranssub))

trans_med2 <- transsizes %>%
  group_by(Focalgene) %>%
  summarize(medeudsub = median(Eudisttranssub))


#Putting back into the original table
for (i in 1:nrow(transsizes)) {
  transsizes[i,"trans_eud_med_sub"] <- trans_med2[trans_med2$Focalgene == transsizes[i,"Focalgene"],"medeudsub"]
  transsizes[i,"trans_eud_avg_sub"] <- trans_avg2[trans_avg2$Focalgene == transsizes[i,"Focalgene"],"avgeudsub"]
}

#plotting the median values for the signature removed
ggplot(data = transsizes, aes(x = Eudistcissub, y = trans_eud_med_sub)) +
  geom_point(color = "black") +
  geom_abline(slope = 1, intercept = 0, color = "seagreen4", size = 2) +
  THEMEMAIN() +
  xlab("Cis-regulatory pleiotropy\n(Euclidean distance)\nSlow growth removed") +
  ylab("Median trans-regulatory pleiotropy\n(Euclidean distance)\nSlow growth removed")
ggsave("Scatterplotcisvstranseudistsubmed.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)

#And now asking what proportion of the focal genes have larger cis-regulatory pleiotropy than trans-reg
length(unique(transsizes[transsizes$trans_eud_med_sub >= transsizes$Eudistcissub,"Focalgene"])) #719 of 748, 96%
nrow(transsizes[transsizes$Eudisttranssub > transsizes$Eudistcissub,]) #6705 of 7162, 94%

#Reading out the final transsizes table
write.table(transsizes, paste0(outputdir,"/transsizeslim.txt"), sep = "\t", row.names = FALSE)

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


#Looking at effect sizes for all focal genes in cis and in trans
for (i in 1:nrow(transsizes)) {
  transsizes[i,"Effectsize"] <- Mvalues[rownames(Mvalues) == transsizes[i,"Strain"],tolower(colnames(Mvalues)) == transsizes[i,"Focalgene"]]
}
for (i in 1:nrow(transsizes)) {
  transsizes[i,"Effectsizecis"] <- Mvalues[rownames(Mvalues) == transsizes[i,"Focalgene"],tolower(colnames(Mvalues)) == transsizes[i,"Focalgene"]]
}

ggplot(data = transsizes, aes(x = Focalgene, y = Effectsize)) +
  geom_point(color = "#336879", alpha = 0.7) +
  geom_point(aes(x = Focalgene, y = Effectsizecis), color = "#f47a60") +
  theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_blank(), legend.position = "none",axis.title = element_text(size = 25), plot.margin = unit(c(2,1,1,1),"cm")) +
  ylab("Focal Gene Expression\n(Log2 Fold Change)") +
  xlab("Focal Gene")
ggsave("Cisvstranseffectall.pdf", plot = last_plot(), path = figdir, width = 10, height = 7)

#Calculating the difference in effect size on the focal gene between cis and trans for each pair to make a histogram
for (i in 1:nrow(transsizes)) {
  transsizes[i,"Effectdiff"] <- 2^(transsizes[i,"Effectsize"]) - 2^(transsizes[i,"Effectsizecis"])
}

ggplot(data = transsizes, aes(x = Effectdiff)) +
  geom_histogram(bins = 100) +
  xlab("Difference in effect size on focal gene\n (trans - cis)") +
  ylab("Count") +
  #xlim(-10,10) +
  THEMEMAIN()
ggsave("Cisvstranseffecthist.pdf", plot = last_plot(), path = figdir, width = 7, height = 7)

ggplot(data = transsizes, aes(x = Effectdiff)) +
  geom_histogram(bins = 100) +
  xlab("Difference in effect size on focal gene\n (trans - cis)") +
  ylab("Count") +
  xlim(-1,10) +
  THEMEMAIN()
ggsave("Cisvstranseffecthistinset.pdf", plot = last_plot(), path = figdir, width = 7, height = 7)

ggplot(data = transsizes, aes(x = Effectsizecis)) +
  geom_histogram(binwidth = 0.5) +
  xlab("Cis-regulatory effects\n on focal gene expression\n(Log2 Fold Change)") +
  ylab("Count") +
  THEMEMAIN()
ggsave("Ciseffecthist.pdf", plot = last_plot(), path = figdir, width = 7, height = 7)

ggplot(data = transsizes, aes(x = Effectsize)) +
  geom_histogram(binwidth = 0.5) +
  xlab("Trans-regulatory effects\n on focal gene expression\n(Log2 Fold Change)") +
  ylab("Count") +
  THEMEMAIN()
ggsave("Transeffecthist.pdf", plot = last_plot(), path = figdir, width = 7, height = 7)
