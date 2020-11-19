if (!requireNamespace("BiocManager", quietly = TRUE))install.packages("BiocManager")
BiocManager::install("GEOquery")
BiocManager::install("multiClust")
BiocManager::install("hgu133plus2a.db")
BiocManager::install("illuminaHumanv4.db")
BiocManager::install("annotate")
BiocManager::install("virtualArray")
BiocManager::install("inSilicoMerging")
# if (!requireNamespace("biocLite.R", quietly = TRUE))install.packages("biocLite.R")
# source("http://bioconductor.org/biocLite.R")
#biocLite("inSilicoMerging")

library(GEOquery)
library(Biobase)
library(multiClust)
library(hgu133plus2a.db)## for othe probe ids
library(illuminaHumanv4.db)### for illumina probe ids
library(annotate)
library(limma)
# library(inSilicoDb)
# library(inSilicoMerging)

# Obtain GSE series matrix file from GEO website using getGEO function
gse <- getGEO(GEO="GSE7410") ### submit GSE id here

# Save the gene expression matrix as an object
data.gse <- exprs(gse[[1]])

# Save the patient clinical data to an object
pheno <- pData(phenoData(gse[[1]]))

# Write the gene expression and clinical data to text files
WriteMatrixToFile(tmpMatrix=data.gse, tmpFileName="/Users/dilrajkaur/Desktop/PANCAN_NETWORK/CESC/GSE7410.expression.txt",
                  blnRowNames=TRUE, blnColNames=TRUE)

WriteMatrixToFile(tmpMatrix=pheno, tmpFileName="/Users/dilrajkaur/Desktop/PANCAN_NETWORK/CESC/GSE7410.clinical.txt",
                  blnRowNames=TRUE, blnColNames=TRUE)


probeset.list <- read.table("/Users/dilrajkaur/Desktop/PANCAN_NETWORK/CESC/GSE7410.expression.txt")
gene.symbols <- getSYMBOL(rownames(probeset.list),"hgu133plus2.db") ### for other type of probe ids 
# gene.symbols <- data.frame(Gene=unlist(mget(x = rownames(probeset.list),envir = illuminaHumanv4SYMBOL,ifnotfound = NA)))###and illuninaHumanv4.db for illumina probe ids--eg.- GSE6253
results <- cbind(gene.symbols,probeset.list)

print(head(results))

############### save results in csv files
write.table(results,"E:/LUNG_comments/GEO_data/GSE6253.expression.csv",row.names=F,col.names=TRUE,sep = ',',append = FALSE)

Z<-read.table(file.choose("E:/LUNG_comments/GEO_data/GSE6253.clinical.txt"), sep="\t")
write.table(Z,"E:/LUNG_comments/GEO_data/GSE6253.clinical.csv",row.names=F,col.names=TRUE,sep = ',',append = FALSE)



