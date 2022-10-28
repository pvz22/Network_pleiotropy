Avalues <- mutsvswt[mutsvswt[,1] == "A",]
Avalues <- Avalues[,2:ncol(Avalues)]
Avalues <- apply(Avalues, c(1,2), as.numeric)
Avalues <- as.data.frame(Avalues, stringsAsFactors = FALSE)
rownames(Avalues) <- mat.df$V1
colnames(Avalues) <- mutantsexwtvar$geneSymbol[2:6124]
colnames(Avalues) <- tolower(colnames(Avalues))

WTvalues <- Avalues/2^(Mvalues)

transsizes <- read.table(paste0(outputdir,"/transsizeslim.txt"), header = 1, sep = "\t")
for (i in 1:nrow(transsizes)) {
  transsizes[i,"Effectsizecis"] <- Mvalues[rownames(Mvalues) == transsizes[i,"Focalgene"],tolower(colnames(Mvalues)) == transsizes[i,"Focalgene"]]
}

for (i in 1:nrow(transsizes)) {
  transsizes[i,"WTmean"] <- mean(WTvalues[,colnames(WTvalues) == transsizes[i,"Focalgene"]], na.rm = TRUE)
}

ggplot(data = transsizes, aes(x = Effectsizecis, y = WTmean)) +
  geom_point() +
  xlab("Deletion log2(fold change)") +
  ylab("Wild Type \nExpression level") +
  THEMEMAIN()

crosshybs <- rownames(transsizes[transsizes$WTmean > 12 & transsizes$Effectsizecis > -2,])
crosshybs <- as.numeric(crosshybs)
transsizesnopara <- transsizes[-crosshybs,]
