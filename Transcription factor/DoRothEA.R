# load required packages
require(viper)
require(data.table)
require(ggplot2)
require(ggpubr)
require(colorspace)

# load transcription factor (TF) regulon gene sets in VIPER format
load('/home/javier/Escritorio/CVID_scBS/Robjects/B_viperRegulon.rdata')

# clean TF names & explore object
names(viper_regulon) = sapply(strsplit(names(viper_regulon), split = ' - '), head, 1)

# explore the regulons object
names(viper_regulon)[1:10]
viper_regulon[[1]]

# load differential expression signature
file <- fread('/home/javier/Escritorio/Covid_autoimmune/DEGs/Noninfected Autoimmune vs Healthy/DEGs_table_MS_Noninfected.csv') 

# exclude probes with unknown or duplicated gene symbols
DEsignature = subset(file, Gene != "" )
DEsignature = subset(DEsignature, ! duplicated(Gene))

# estimate z-score values for the GES (check VIPER manual for details)
myStatistics = matrix(DEsignature$logFC, dimnames = list(DEsignature$Gene, 'logFC') )
myPvalue = matrix(DEsignature$adj.P.Val, dimnames = list(DEsignature$Gene, 'P.Value') )
mySignature = (qnorm(myPvalue/2, lower.tail = FALSE) * sign(myStatistics))[, 1]
mySignature = mySignature[order(mySignature, decreasing = T)]

# estimate TF activities
mrs = msviper(ges = mySignature, regulon = viper_regulon, minsize = 4, ges.filter = F)
TF_activities = data.frame(Regulon = names(mrs$es$nes),
                           Size = mrs$es$size[ names(mrs$es$nes) ],
                           NES = mrs$es$nes,
                           p.value = mrs$es$p.value,
                           FDR = p.adjust(mrs$es$p.value, method = 'fdr'))
TF_activities = TF_activities[ order(TF_activities$p.value), ]

# plot results
plot = TF_activities[TF_activities$FDR < 0.05,]
plot = plot[order(plot$NES, decreasing = T),]
plot$Regulon = reorder(plot$Regulon, plot$NES*-1)
ggplot(plot, aes(fill=-log10(FDR), x=Regulon, y=NES, color="black")) + 
  geom_col(color="black", size=0.2) + theme_bw() + scale_fill_gradient(high = "#132B43", low = "#56B1F7") + coord_flip()