#AdrienneVancura June 2018
#prepare barplot for showing hCIS distribution among gene families
#install.packages("forcats")
#run script within bash script
library("ggplot2")
library("RColorBrewer")
library("plyr")
library("dplyr")
library("forcats")

#https://stackoverflow.com/questions/5247766/running-an-r-file-with-arguments-from-a-bash-script
#args <- commandArgs()
#args = commandArgs(trailingOnly=TRUE)


df <- read.table ("CISelementsGeneFamiliesPercentages.txt", header=F)
df$specie <- c("CLC","CLC","CLC","CLC","CLC")

hCIS <- df$V2
geneclass <- df$V1
CLC <- df$specie

#bp <- ggplot(df, aes(1, hCIS, fill=geneclass))+ geom_bar(stat="identity") 
#bp + scale_fill_manual(values=c('#a8ddb5','#7bccc4', '#6baed6','#fb6a4a','#a50f15'))
#ggsave("GeneTypes.pdf", width=5, height=7)


#
#
#
#
#


df %>%
  mutate(geneclass = fct_relevel(geneclass, "CLC2", "nonCLC", "CGC", "nonCGC", "intergenic")) %>%
  ggplot( aes(x=CLC, y=hCIS, fill=geneclass)) +
  scale_fill_manual(values=c('#7bccc4','#bae4bc','#0868ac','#43a2ca','#cbc9e2')) +
  geom_bar(stat="identity") + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ggsave("GeneTypes.pdf", width=5, height=7)












