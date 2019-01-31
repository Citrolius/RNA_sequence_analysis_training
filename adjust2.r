library(DESeq2)

setwd("//salt/bioinf_training/colin/RNA-seq/sort_bam_files/")
getwd()

counts <-read.table("T:/colin/RNA-seq/sort_bam_files/counts.txt",row.names=1, header=TRUE, sep='\t')
head(counts) ; dim(counts)

#design1 <- c(rep("R",4), rep("S",4))
#design1 <- c(rep("R",4), rep("S",3))
#design1 <- c(rep("R",4), rep("S",2))
design1 <- c(rep("R",3), rep("S",2))

colData <- data.frame(design1)
colData

#dds <- DESeqDataSetFromMatrix(countData=counts, design=~design1,colData=colData)
#dds <- DESeqDataSetFromMatrix(countData=countsr6, design=~design1,colData=colData)
#dds <- DESeqDataSetFromMatrix(countData=countsr7, design=~design1,colData=colData)
dds <- DESeqDataSetFromMatrix(countData=countsr8, design=~design1,colData=colData)

dds <- DESeq(dds)
dds

res <- results(dds)
res ; summary(res)

#res.ordered <- res[order(res$padj),]
res.ordered <- res[order(res$pvalue),]
res.ordered ; summary(res.ordered)

tmpdir <-"tmpdir_temporary"
tmpath <-file.path("//salt/bioinf_training/colin/RNA-seq/sort_bam_files/")
dir.create(tmpath, showWarnings=FALSE)


file <-paste(tmpath,"differential_expression3.csv",sep='/')
write.csv(as.data.frame(res.ordered),file=file)

#log transform and plot PCA
rld <- rlog(dds, blind=FALSE)
rld
plotPCA(rld, intgroup="design1") 

#Get coordinates of PCA plots
plotPCA(rld, intgroup="design1", returnData=TRUE) 
###################################################################################################################################################

#NW6 removed
countsr6 <- counts[,-6]

#NW8 removed
countsr7 <- countsr6[,-7]

#NW4 removed
countsr8 <- countsr7[,-4]





