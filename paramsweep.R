library(dplyr)
library(cowplot)
library(ggplot2)
outputdir <- "/Users/petravandezande/Documents/Output/Projects/Pleiotropy/Friend"
figdir <- "/Users/petravandezande/Documents/Figures/Projects/Pleiotropy/Friend"

Mvalues <- read.table(paste0(outputdir,"/Mvalues.txt"), sep="\t")
pvalues <- read.table(paste0(outputdir,"/p.values.txt"), sep = "\t")

Pleiotropy <- matrix(NA, nrow = nrow(Mvalues), ncol = ncol(Mvalues))
a <- c(0.2630344,0.5849625, 1, 1.584963) #Possible Mvalues
b <- c(0.2,0.1,0.05,0.01,0.001) #Possible p-values

for (i in 1:length(a)) {
  for (j in 1:length(b)) {
    X <- matrix(NA, nrow = nrow(Mvalues), ncol = ncol(Mvalues))
    X <- ifelse(abs(Mvalues) >= a[i] & pvalues <= b[j], 1, 0)
    X <- assign(paste0("Pleiotropy",i,j), X)
  }
}

for (k in c(1:4)) {
  for (l in c(1:5)) {
    Pleiotropy <- get(paste0("Pleiotropy",k,l))
    colnames(Pleiotropy) <- colnames(Mvalues)
    rownames(Pleiotropy) <- rownames(Mvalues)
    colnames(Pleiotropy) <- tolower(colnames(Pleiotropy))
    cissize <- c()
    transsizes <- data.frame()
    for (i in 1:nrow(Pleiotropy)) {
      if (rownames(Pleiotropy)[i] %in% colnames(Pleiotropy)) {
        print(i)
        cissize[i] <- rowSums(Pleiotropy)[i]
        focalgene <- rownames(Pleiotropy)[i]
        if(sum(Pleiotropy[,focalgene]) > 1) { #To eliminate where it is only the self loop
          transregs <- as.data.frame(Pleiotropy[Pleiotropy[,focalgene] == 1,])
          transregs$Size <- rowSums(transregs)
          trans <- data.frame("Strain" = rownames(transregs),"Size" = transregs$Size, "Focalgene" = focalgene)
          transsizes <- rbind(transsizes, trans)
        }
      }
    }
    cissize <- assign(paste0("cissize",k,l), cissize)
    transsizes <- assign(paste0("transsizes",k,l), transsizes)
  }
}


for (k in c(1:4)) {
  for (l in c(1:5)) {
    cissize <- get(paste0("cissize",k,l))
    transsizes <- get(paste0("transsizes",k,l))
    cissize <- cissize[1:1484] #Getting rid of these other environments
    cissize <- data.frame("Cissize" = cissize)
    rownames(cissize) <- rownames(Pleiotropy)[1:1484]
    transsizes$Focalgene <- as.character(transsizes$Focalgene)
    transsizes$Strain <- as.character(transsizes$Strain)
    #Removing self-loops as trans-regulatory perturbations
    transsizes <- transsizes[which(transsizes$Strain != transsizes$Focalgene),]
    trans_avg <- transsizes %>%
      group_by(Focalgene) %>%
      summarize(avg = mean(Size))
    
    trans_med <- transsizes %>%
      group_by(Focalgene) %>%
      summarize(med = median(Size))
    
    trans_avg <- as.data.frame(trans_avg)
    for (i in 1:nrow(trans_avg)) {
      trans_avg[i,"Difference"] <- trans_avg[i,"avg"] - cissize[as.character(trans_avg[i,"Focalgene"]),"Cissize"]
    }
    trans_med <- as.data.frame(trans_med)
    for (i in 1:nrow(trans_med)) {
      trans_med[i,"Difference"] <- trans_med[i,"med"] - cissize[as.character(trans_med[i,"Focalgene"]),"Cissize"]
    }
    trans_avg <- assign(paste0("trans_avg",k,l), trans_avg)
    trans_avg <- assign(paste0("trans_med",k,l), trans_med)
  }
}

n <- 1
for (k in 1:4) {
  for (l in 1:5) {
    trans_avg <- get(paste0("trans_avg",k,l))
    plot <- ggplot(data = trans_avg, aes(x=Difference)) +
      geom_histogram() +
      geom_vline(xintercept = 1)
    plot <- assign(paste0("p",n), plot)
    n <- n+1
  }
}
plotlist <- list(p1,  p2,  p3,  p4,  p5,  p6,  p7,  p8,  p9,  p10, p11, p12, p13,
              p14, p15, p16, p17, p18, p19, p20)
plot_grid(plotlist = plotlist, nrow = 4, ncol = 5)
ggsave("paramsweep.pdf",plot = last_plot(), path = figdir, width = 12, height = 12)


#Limiting the subsetting stuff to only genes that are reducing expression of the target gene
Pleiotropydown <- matrix(NA, nrow = nrow(Mvalues), ncol = ncol(Mvalues))
Pleiotropydown <- ifelse(Mvalues <= -0.765  & pvalues <= 0.05, 1, 0)
colnames(Pleiotropydown) <- colnames(Mvalues)
rownames(Pleiotropydown) <- rownames(Mvalues)
colnames(Pleiotropydown) <- tolower(colnames(Pleiotropydown))
cissize <- c()
transsizes <- data.frame()
for (i in 1:nrow(Pleiotropydown)) {
  if (rownames(Pleiotropydown)[i] %in% colnames(Pleiotropydown)) {
    print(i)
    cissize[i] <- rowSums(Pleiotropydown)[i]
    focalgene <- rownames(Pleiotropydown)[i]
    if (sum(Pleiotropydown[,focalgene]) > 1) { #To eliminate where it is only the self loop
      transregs <- as.data.frame(Pleiotropydown[Pleiotropydown[,focalgene] == 1,])
      transregs$Size <- rowSums(transregs)
      trans <- data.frame("Strain" = rownames(transregs),"Size" = transregs$Size, "Focalgene" = focalgene)
      transsizes <- rbind(transsizes, trans)
    }
  }
}
#Hopped over to the "deletiondata.R" script to take a look at how this compared for the other analyses. Same patterns


for (i in 1:nrow(transsizes)) {
  if (transsizes[i,"cissize"] > 0) {
    fgvec <- which(Pleiotropydown[transsizes[i,"Focalgene"],] == 1)
    transsizes[i,"fracnest"] <- sum(Pleiotropydown[transsizes[i,"Strain"],fgvec])/length(fgvec) 
  }
}
write.table(transsizes, paste0(outputdir,"/transsizes.txt"), sep = "\t")

table(transsizes$fracnest >= 1)
table(transsizes[transsizes$fracnest >= 1,"cissize"] == 1)
table(transsizes[transsizes$cissize > 1,"fracnest"] == 1) #118 true, 1318 false
table(transsizes[transsizes$cissize == 1,"fracnest"] == 1) #6 false, 1217 true

table(table(transsizes[transsizes$cissize != 1 & transsizes$fracnest == 1,"Strain"]))
table(transsizes[transsizes$cissize != 1 & transsizes$fracnest == 1,"Strain"])

pie(c(1223, 118, 1318), labels = c("Cis=1","Prop=1","Prop<1"), col = c("white", "lightgrey","darkgrey"))

ggplot(data = transsizes[transsizes$cissize != 1,], aes(x = fracnest)) + 
  geom_histogram(fill = "darkgrey", bins = 10) +
  THEMEMAIN() +
  xlab("Proportion of cis-regulatory DEGs\npresent in trans-regulatory DEGs") +
  ylab("Count")
ggsave("Nestpropdown.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)



