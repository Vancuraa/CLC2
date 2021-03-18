install.packages('RMySQL', repos='http://cran.us.r-project.org')

local({r <- getOption("repos")
r["CRAN"] <- "http://cran.r-project.org" 
options(repos=r)
})

chooseCRANmirror(graphics=FALSE, ind=55)

data=read.table("OUTPUTforBarplot.txt",header=T)

install.packages("edgeR")
library(edgeR) 

install.packages("ggplot2")
library(ggplot2) 
library(reshape2)

source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
biocLite("forcats") 


ggplot(data=data, aes(fill=Gene_type, y=Value,x=Feature))+
  geom_bar(stat="identity", position="fill")+
  #geom_text(data=data,aes(label=data$hCIS_elements), position="fill", vjust=-0.25)+
  #labs(x = "Gene type", y = "hCIS / Mb")+
  theme_bw()+
  scale_fill_manual(values=c("deepskyblue1","deepskyblue3","firebrick1","firebrick3","gray75")) +
  theme(panel.grid=element_blank())
ggsave("plot1.pdf", width = 6, height = 4)


allNucleotides <- data[1:5, ]
hCISelements <- data[6:10, ]
total <- merge(allNucleotides,hCISelements,by="Gene_type")
Calculated <- transform(total, new = Value.y / Value.x * 1000000)
table <- data.frame(Calculated$Gene_type,Calculated$Value.x,Calculated$Value.y, Calculated$new)
colnames(table) <- c("Gene_type", "All_Nucleotides", "hCIS_elements", "hCIS_per_nucleotide")
data<-table


#data=read.table("data_for_R_barplot2.txt",header=T)
ggplot(data=data, aes(fill=Gene_type, y=hCIS_per_nucleotide,x=Gene_type))+
  geom_bar(stat="identity")+
  geom_text(data=data,aes(label=data$hCIS_elements), vjust=-0.25)+
  #labs(x = "Gene type", y = "hCIS / Mb")+
  theme_bw()+
  scale_fill_manual(values=c("deepskyblue1","deepskyblue3","firebrick1","firebrick3","gray75")) +
  theme(panel.grid=element_blank())
ggsave("plot2.pdf", width = 6, height = 4)




