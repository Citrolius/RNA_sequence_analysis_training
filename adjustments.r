# Analysing RNA-sequencing results using differential expression with DESEq2
#
###########################################################################

library(DESeq2)
#DESeq2 estimates mean dependence in count data.
#It uses negative binomial distribution to test differential expression.


setwd("//salt/bioinf_training/colin/RNA-seq/sort_bam_files/")
getwd()

#########################################################################################################

#1) Load the count table you generated with featureCounts then show the heads and dimensions of the data.
counts <-read.table("T:/colin/RNA-seq/sort_bam_files/counts.txt",row.names=1, header=TRUE, sep='\t')
head(counts) ; dim(counts) 
#design1 <- c(rep("R",4),rep("S",4))
#design1 <- c(rep("R",4), rep("S",3))
#design1 <- c(rep("R",3), rep("S",3))
design1 <- c(rep("R",3), rep("S",2))



#2) Create a small dataframe in order to teach DESEq2 how to interpret the columns.
#You need to tell DESeq2 that column data is the previous variable you have set.
#You already told it the row name in step 1. These are the formats DESeq2 must know.
colData <- data.frame(design1)
colData


#3) Create DESeqDataSet to format the data so DESeq2 can read it.
# DESeqDataSetFromMatrix MUST HAVE A DESIGN TO WORK!!!!!!!
#dds <- DESeqDataSetFromMatrix(countData=counts, design=~design1,colData=colData)
#dds <- DESeqDataSetFromMatrix(countData=counts_no_outliers, design=~design1,colData=colData)
#dds <- DESeqDataSetFromMatrix(countData=counts_no_outliers2, design=~design1,colData=colData)
dds <- DESeqDataSetFromMatrix(countData=counts_no_outliers3, design=~design1,colData=colData)
dds

#4) Basic DESEq2 analysis
dds <- DESeq(dds)
dds

#5) Extract and sort results by p-values.

res <- results(dds)
res ; summary(res)

#res.ordered <- res[order(res$padj),]
res.ordered <- res[order(res$pvalue),]
res.ordered ; summary(res.ordered)

#6) export to temporary directory
tmpdir <-"tmpdir_temporary"
tmpath <-file.path("//salt/bioinf_training/colin/RNA-seq/sort_bam_files/")
dir.create(tmpath, showWarnings=FALSE)

file <-paste(tmpath,"differential_expression.csv",sep='/')
write.csv(as.data.frame(res.ordered),file=file)

###########################################################################
#Use a MA plot to compare Resistance and succeptible perhaps?
#plotMA( res, ylim = c(-1, 1) )


rld <- rlog(dds, blind=FALSE)
rld
plotPCA(rld, intgroup="design1") 
plotPCA(rld, intgroup="design1", returnData=TRUE) 
#MW5 removed
counts_no_outliers <- counts[,-5]
#MW4 removed
counts_no_outliers2 <- counts_no_outliers[,-4]
#MW8 removed
counts_no_outliers3 <- counts_no_outliers2[,-6]



#MW6 removed instead of MW8
counts_no_outliers3 <- counts_no_outliers2[,-4]


