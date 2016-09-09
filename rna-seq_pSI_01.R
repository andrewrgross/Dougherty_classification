### Tissues Specific Expression Analysis, Andrew R Gross, 2016-08-26
### Using a package developed by the Dougherty lab, calculate the uniqueness of a set of genes and the probability
### of another set containing those genes.

########################################################################
### Header
########################################################################
library(biomaRt)
library(reshape2)
library(gdata)
library(ggplot2)
library(pSI)
load("C:/Users/grossar/Bioinform/DATA/rna_seq/Dougherty/pSI.data/data/human.rda")
help(pSI)

ensembl = useMart(host='www.ensembl.org',biomart='ENSEMBL_MART_ENSEMBL',dataset="hsapiens_gene_ensembl")
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)

########################################################################
### Functions
########################################################################

addMedSD <- function(dataframe) {
  median <- apply(dataframe,1,median)
  sd <- apply(dataframe,1,sd)
  return(data.frame(dataframe,median,sd))  
}
sortByMed <- function(dataframe) {
  order <- order(dataframe$median,decreasing=TRUE)
  return(dataframe[order,])
}
convertIDs <- function(dataframe) {
  ensemblIDs <- c()
  for (rowName in row.names(dataframe)) {
    ensemblID <- strsplit(rowName,"\\.")[[1]][1]
    ensemblIDs <- c(ensemblIDs,ensemblID)
  }
  row.names(dataframe)<-ensemblIDs
  return(dataframe)
}
addGene <- function(dataframe) {
  genes <- getBM(attributes=c('ensembl_gene_id','external_gene_name'), filters='ensembl_gene_id', values=row.names(dataframe), mart=ensembl)
  genes <- genes[match(row.names(dataframe),genes[,1]),]
  Gene <- c()
  for (rowNumber in 1:length(genes[,1])) {
    newGene <- genes[rowNumber,][,2]
    Gene <- c(Gene, newGene)
  }
  dataframe[length(dataframe)+1] <- Gene
  names(dataframe)[ncol(dataframe)] <- "Gene"
  return(dataframe)
}
reorder_size <- function(x) {
  factor(x, levels = names(sort(table(x))))
}
########################################################################
### Data input
########################################################################

# Normalized refseq data in units of TPM
TPMdata <- read.csv("z://Data/RNAseq HT neurons and tissue/Andrews_files/20160317_Sareen_rerun-29299281.tpm.csv", row.names=1)

# Metadata
metadata.sam <- read.csv("z://Data/RNAseq HT neurons and tissue/2nd rerun/samples.csv", row.names=1)

### Gtex references
references <- read.table("c://Users/grossar/Bioinform/DATA/rna_seq/Reference_transcriptomes/GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct",sep="\t",header=TRUE,row.names=1)

########################################################################
### Formatting
########################################################################
### Format TPM
TPMdata <- TPMdata[row.names(metadata.sam)]  # Reorder columns to match metadata

### Remove unwanted column
TPMdata <- TPMdata[-match("Reference",names(TPMdata))]
metadata.sam <- metadata.sam[-match("Reference",row.names(metadata.sam)),]

########################################################################
### pSI determination
########################################################################

#pSI.input <- sample.data$pSI.input
pSI.input <- references[2:ncol(references)]
pSI.output <- specificity.index(pSI.input)
pSI.thresholds <- pSI.list(pSI.output)

########################################################################
### Save pSI tables
########################################################################

setwd("z://Data/RNAseq HT neurons and tissue/Andrews_files/")
write.csv(pSI.output,"pSI_gTEX_full_8-26-16.csv")

########################################################################
### Load pSI tables
########################################################################

### Dougherty tissue list, NAR 2015, Table S3 -- 25 tissues, 18056 genes labeled with gene names
psi.tsea<-read.table("http://genetics.wustl.edu/jdlab/files/2015/10/TableS3_NAR_Dougherty_Tissue_gene_pSI_v3.txt",header=T,row.names=1)

### Dougherty cell types, JoN 2016 -- 60 cell types, 16,866 genes in gene name format
psi.csea <- human$developmental.periods$psi.out

### Gtex tissues -- 53 tissues, 56,418 genes labeled with ENSB IDs
psi.gtex <- read.csv("z:/Data/RNAseq HT neurons and tissue/Andrews_files/pSI_gTEX_full_8-26-16.csv", header=T,row.names=1)
psi.gtex <- convertIDs(psi.gtex)

########################################################################
### Subsetting
########################################################################

metadata.sam.iHT <- metadata.sam[metadata.sam$Group == 'iHT',]  # Select metadata from iHT samples
metadata.sam.iHT.f <- metadata.sam.iHT[metadata.sam.iHT$Sex =='F',]  # Select female iHT samples
metadata.sam.HT <- na.omit(metadata.sam[metadata.sam$Type == "HT",])  # Select HT samples
metadata.sam.HT.f <- metadata.sam.HT[metadata.sam.HT$Sex =='F',]  # Select female iHT samples
metadata.sam.iPSC <- metadata.sam[metadata.sam$Source == "iPSC",]  # Select iPSC samples
metadata.sam.aHT <- metadata.sam[metadata.sam$Source == "Adult",]  # Select adult HT samples
metadata.sam.MN <- na.omit(metadata.sam[metadata.sam$Type == "MN",])  # Select HT samples

HT.data <- TPMdata[match(row.names(metadata.sam.HT),names(TPMdata))]
iHT.data <- TPMdata[match(row.names(metadata.sam.iHT),names(TPMdata))]
aHT.data <- TPMdata[match(row.names(metadata.sam.aHT),names(TPMdata))]
iPSC.data <- TPMdata[match(row.names(metadata.sam.iPSC),names(TPMdata))]
MN.data <- TPMdata[match(row.names(metadata.sam.MN),names(TPMdata))]
iHT.f.data <- TPMdata[match(row.names(metadata.sam.iHT.f),names(TPMdata))]
HT.f.data <- TPMdata[match(row.names(metadata.sam.HT.f),names(TPMdata))]

datalist <- list(HT.data,iHT.data,aHT.data,iPSC.data,MN.data,iHT.f.data,HT.f.data)
names(datalist) <- c("HT","iHT","aHT","iPSC","MN","iHT.f","HT.f")

########################################################################
### Reorder genes by expression
########################################################################
### Add Median and sort by median
for(list.element in 1:length(datalist)) {
  datalist[[list.element]] <- sortByMed(addMedSD(datalist[[list.element]]))[1:2000,]
}

### Generate full lists of genes
genelist <- datalist
for(list.element in 1:length(datalist)) {
  current.genelist <- addGene(genelist[[list.element]])
  genelist[[list.element]] <- current.genelist[ncol(current.genelist)]
}

########################################################################
### Subset gene list
########################################################################

names(genelist)
geneset <- 2
genes.of.interest.full <- genelist[[geneset]][,1]
ids.of.interest.full <- row.names(genelist[[geneset]])
geneset <- names(genelist[geneset])
print(geneset)

length <- 2200

genes.of.interest <- genes.of.interest.full[1:length]
ids.of.interest <- ids.of.interest.full[1:length]


########################################################################
### Calculate overlap p-values
########################################################################

tsea.results <- fisher.iteration(pSIs = psi.tsea, candidate.genes = genes.of.interest)
csea.results <- fisher.iteration(pSIs = psi.csea, candidate.genes = genes.of.interest)
gtex.results <- fisher.iteration(pSIs = psi.gtex, candidate.genes = ids.of.interest)

########################################################################
### Report significant genes
########################################################################


########################################################################
### Save p-value results
########################################################################

timestamp <- strftime(Sys.time(),"%a%b%d%H%M")
setwd("z://Data/RNAseq HT neurons and tissue/Andrews_files/pSI_files/pSI_results/")

write.csv(tsea.results, paste0("tsea_",geneset,"_",length,"_",timestamp,".csv"))
write.csv(csea.results, paste0("csea_",geneset,"_",length,"_",timestamp,".csv"))
write.csv(gtex.results, paste0("gtex_",geneset,"_",length,"_",timestamp,".csv"))

########################################################################
### Load p-value results
########################################################################

#setwd("z://Data/RNAseq HT neurons and tissue/Andrews_files/pSI_files/pSI_results/")

#tsea.results <- read.csv("tsea_")
#csea.results <- read.csv("csea_")
#gtex.results <- read.csv("gtex_")

########################################################################
### Visualize results
########################################################################
### Plot bar graphs 

### Format TSEA data
tsea.plotting <- tsea.results[tsea.results$`0.05 - adjusted` != 1,]
names(tsea.plotting) <- c('p0.05','p0.01','p0.001','p1e-4')
tsea.plotting <- -log10(tsea.plotting)
tsea.plotting$tissue <- row.names(tsea.plotting)
tsea.plotting.melt <- melt(tsea.plotting,id.var = "tissue")

### Plot TSEA data
tsea.plot <- ggplot(data = tsea.plotting.melt,aes(x = reorder(tissue, -value), y = value, fill = variable)) +
  geom_bar(position = "dodge",stat="identity") + 
  labs(title = paste("p-value of the overlap of the top",length,"genes in",geneset),
       x = "Tissue",
       y = "Negative log10 p-value") +
  theme(plot.title = element_text(color="black", face="bold", size=22, margin=margin(10,0,20,0)),
        axis.title.x = element_text(face="bold", size=14,margin =margin(20,0,10,0)),
        axis.title.y = element_text(face="bold", size=14,margin =margin(0,20,0,10)),
        plot.margin = unit(c(1,1,1,1), "cm"), axis.text = element_text(size = 12),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_brewer("p-value",palette = 'OrRd', direction = -1)

tsea.plot

### Format CSEA data
csea.plotting <- csea.results[csea.results$`0.05 - adjusted` != 1,]
names(csea.plotting) <- c('p0.05','p0.01','p0.001','p1e-4')
csea.plotting <- -log10(csea.plotting)
csea.plotting$tissue <- row.names(csea.plotting)
csea.plotting.melt <- melt(csea.plotting,id.var = "tissue")

### Plot CSEA data
csea.plot <- ggplot(data = csea.plotting.melt,aes(x = reorder(tissue, -value), y = value, fill = variable)) +
  geom_bar(position = "dodge",stat="identity") + 
  labs(title = paste("p-value of the overlap of the top",length,"genes in",geneset),
       x = "Tissue",
       y = "Negative log10 p-value") +
  theme(plot.title = element_text(color="black", face="bold", size=22, margin=margin(10,0,20,0)),
        axis.title.x = element_text(face="bold", size=14,margin =margin(20,0,10,0)),
        axis.title.y = element_text(face="bold", size=14,margin =margin(0,20,0,10)),
        plot.margin = unit(c(1,1,1,1), "cm"), axis.text = element_text(size = 12),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_brewer("p-value",palette = 'OrRd', direction = -1)

csea.plot

### Format gtex data
gtex.plotting <- gtex.results[gtex.results$`0.01 - adjusted` != 1,]
names(gtex.plotting) <- c('p0.05','p0.01','p0.001','p1e-4')
gtex.plotting <- gtex.plotting[c('p0.01','p0.001','p1e-4')]
gtex.plotting <- -log10(gtex.plotting)
gtex.plotting$tissue <- row.names(gtex.plotting)
gtex.plotting[order(gtex.plotting$p0.01,decreasing = TRUE),]
gtex.plotting.melt <- melt(gtex.plotting,id.var = "tissue")

### Plot gtex data
gtex.plot <- ggplot(data = gtex.plotting.melt,aes(x = reorder(tissue, -value), y = value, fill = variable)) +
  geom_bar(position = "dodge",stat="identity") + 
  labs(title = paste("p-value of the overlap of the top",length,"genes in",geneset),
       x = "Tissue",
       y = "Negative log10 p-value") +
  theme(plot.title = element_text(color="black", face="bold", size=22, margin=margin(10,0,20,0)),
        axis.title.x = element_text(face="bold", size=14,margin =margin(20,0,10,0)),
        axis.title.y = element_text(face="bold", size=14,margin =margin(0,20,0,10)),
        plot.margin = unit(c(1,1,1,1), "cm"), axis.text = element_text(size = 12),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_brewer("p-value",palette = 'Reds')

gtex.plot

########################################################################
### Output plots
########################################################################

timestamp <- strftime(Sys.time(),"%a%b%d%H%M")
setwd("z://Data/RNAseq HT neurons and tissue/Andrews_files/pSI_files/pSI_results/")

tiff(filename=paste0("tsea_",geneset,"_",length,"_",timestamp,".tiff"), 
     type="cairo",
     units="in", 
     width=14, 
     height=14, 
     pointsize=12, 
     res=100)
tsea.plot
dev.off()

tiff(filename=paste0("csea_",geneset,"_",length,"_",timestamp,".tiff"), 
     type="cairo",
     units="in", 
     width=14, 
     height=14, 
     pointsize=12, 
     res=100)
csea.plot
dev.off()


tiff(filename=paste0("gtex_",geneset,"_",length,"_",timestamp,".tiff"), 
     type="cairo",
     units="in", 
     width=14, 
     height=14, 
     pointsize=12, 
     res=100)
gtex.plot
dev.off()





###################################################
### Generate overlap lists
###################################################

test <- candidate.overlap(pSIs=psi_values,candidate.genes=test.genes)

test <- candidate.overlap(brain_genes_0.05,genes.of.interest)



