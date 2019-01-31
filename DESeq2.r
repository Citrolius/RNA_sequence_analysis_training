############################################################################
# SCRIPT TO RUN DEseq2 on HISAT2 output
# ======================================
###########################################################################

#This script analyses RNA-sequencing results using differential expression with DESEq2
#It uses negative binomial distribution to test differential expression.
###########################################################################
Things to note
==============
#1) Load the count table you generated with featureCounts.
#You must state the first row is the row name as Deseq2 doesn't like headers in data.
#you must create a variable that defines the structure of the experiment. i.e design1, the list consists of 2 replicates, R and S.

#2) Create a small dataframe in order to teach DESEq2 how to interpret the columns.
#You need to tell DESeq2 that column data is the previous variable you have set.
#You already told it the row name in step 1. These are the formats DESeq2 must know.

#3) Create DESeqDataSet to format the data so DESeq2 can read it.
# DESeqDataSetFromMatrix MUST HAVE A DESIGN TO WORK!!!!!!!





###########################################################################

#load DESeq2 package
library(DESeq2)

#set working directory
setwd("//salt/bioinf_training/colin/RNA-seq/sort_bam_files/")
getwd()




counts <-read.table("T:/colin/RNA-seq/sort_bam_files/counts.txt",row.names=1, header=TRUE, sep='\t')
head(counts) ; dim(counts) 
design1 <- c(rep("R",4),rep("S",4)) 




colData <- data.frame(design1)
colData


dds <- DESeqDataSetFromMatrix(countData=counts, design=~design1,colData=colData)
#dds <- DESeqDataSetFromMatrix(countData=counts_no_outliers, design=~design1,colData=colData)
dds

dds <- DESeq(dds)
dds

#Extract and sort results by p-values.
res <- results(dds)
res ; summary(res)

#res.ordered <- res[order(res$padj),]
res.ordered <- res[order(res$pvalue),]
res.ordered ; summary(res.ordered)

#export to temporary directory
tmpdir <-"tmpdir_temporary"
tmpath <-file.path("//salt/bioinf_training/colin/RNA-seq/sort_bam_files/")
dir.create(tmpath, showWarnings=FALSE)

file <-paste(tmpath,"differential_expression.csv",sep='/')
write.csv(as.data.frame(res.ordered),file=file)

#log transform the data.
rld <- rlog(dds, blind=FALSE)
rld

#plot Principal Component Analysis.
plotPCA(rld, intgroup="design1") 

#get coordinates
plotPCA(rld, intgroup="design1", returnData=TRUE) 


