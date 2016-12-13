### pSI tests -- Andrew R Gross , 2016-11-28
### A set of scripts intended to better understand the pSI classification system.

########################################################################
### Header
########################################################################
library(pSI)
load("C:/Users/grossar/Bioinform/DATA/rna_seq/Dougherty/pSI.data/data/human.rda")
help(pSI)

ensembl = useMart(host='www.ensembl.org',biomart='ENSEMBL_MART_ENSEMBL',dataset="hsapiens_gene_ensembl")
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)

########################################################################
### Functions
########################################################################

########################################################################
### Import Data
########################################################################
### Dougherty tissue list, NAR 2015, Table S3 -- 25 tissues, 18056 genes labeled with gene names
psi.tsea<-read.table("http://genetics.wustl.edu/jdlab/files/2015/10/TableS3_NAR_Dougherty_Tissue_gene_pSI_v3.txt",header=T,row.names=1)

exprs.human <- human$young.adulthood$psi.in
psi.human <- human$young.adulthood$psi.out

########################################################################
### Format
########################################################################
summary(exprs.human)
exprs.human.log <- log2(exprs.human+1)

### Generating pSI tables to compare to the provided one
psi.test.log <- specificity.index(exprs.human.log)
psi.test <- specificity.index(exprs.human)
psi.test.log.2 <- specificity.index(exprs.human.log,e_min=1.3)

# Check first few columns.  Count the NAs.
head(psi.human)
head(psi.test)
head(psi.test.log)

psi.human[10:30,]
psi.test[10:30,]
psi.test.log[10:30,]
psi.test.log.2[10:30,]

test.list <- list()
test.list <- list(test.list,psi.test)

