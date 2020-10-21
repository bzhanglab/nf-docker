library(corrplot)
library(purrr)
library(tidyr)
library(ggplot2)
library(preprocessCore)
library(edgeR)
library(factoextra)

rm(list=ls())
output_dir <- "RSEM_results_summarized"
dir.create(output_dir)
matrix_name_flag="gene_rsem_with_circRNA"
sample_list_for_summary=list.dirs(path = ".", full.names = FALSE, recursive = FALSE)
sample_list_for_summary=sample_list_for_summary[sample_list_for_summary!=output_dir]
example_file_path=file.path(sample_list_for_summary[1], paste0(matrix_name_flag,"_BCM.txt"))
example_table=read.table(example_file_path,header=TRUE,row.names=1,check.name=FALSE,sep="\t")
total_item_number=nrow(example_table)
raw_matrix=matrix(data=0,nrow=total_item_number,ncol=length(sample_list_for_summary))
normalized_matrix=matrix(data=0,nrow=total_item_number,ncol=length(sample_list_for_summary))
rownames(raw_matrix)=rownames(example_table)
colnames(raw_matrix)=sample_list_for_summary

for(i in 1:length(sample_list_for_summary))
{
    input_file_path=file.path(sample_list_for_summary[i], paste0(matrix_name_flag,"_BCM.txt"))
    print(input_file_path)
    sample_table=read.table(input_file_path,header=TRUE,row.names=1,check.name=FALSE,sep="\t")
    raw_matrix[,i]=sample_table[,1]
}

#normalized_matrix <- normalize.quantiles(as.matrix(raw_matrix), copy = TRUE)
#row.names(normalized_matrix) <- row.names(raw_matrix)
#colnames(normalized_matrix) <- colnames(raw_matrix)

group <- as.factor(colnames(raw_matrix))
cds <- DGEList(counts = raw_matrix, group=group)
cds <- calcNormFactors(cds, method="upperquartile")
scale <- cds$samples$lib.size*cds$samples$norm.factors
normalized_matrix <- round(t(t(raw_matrix)/scale)*mean(scale),2)

write.table(round(log2(raw_matrix+1),2),file=file.path(output_dir, paste0(matrix_name_flag,"_tumor_normal_raw_log2(x+1)_BCM.txt")),quote = FALSE,sep="\t")
write.table(round(log2(normalized_matrix+1),2),file=file.path(output_dir, paste0(matrix_name_flag,"_tumor_normal_UQ_log2(x+1)_BCM.txt")),quote = FALSE,sep="\t")

pdf(file=file.path(output_dir, paste0(matrix_name_flag,"_tumor_normal_raw_log2(x+1)_BCM_boxplot.pdf")),height=7,width=12)
boxplot(log2(raw_matrix+1),las = 2,ylab="Log2(RSEM+1)",main="Distribution of gene expression")
dev.off()
pdf(file=file.path(output_dir, paste0(matrix_name_flag,"_tumor_normal_UQ_log2(x+1)_BCM_boxplot.pdf")),height=7,width=12)
boxplot(log2(normalized_matrix+1),las = 2,ylab="Log2(RSEM+1)",main="Distribution of gene expression")
dev.off()

raw_matrix=raw_matrix[rowMeans(raw_matrix)>0,]
M <- cor(as.matrix(raw_matrix),method="spearman")
pdf(file=file.path(output_dir, paste0(matrix_name_flag,"_tumor_normal_raw_log2(x+1)_BCM_corplot.pdf")),height=7,width=7)
corrplot(M, method = "color",order = "hclust",tl.cex = 0.3,cl.lim=c(min(M),max(M)), is.corr=FALSE)
dev.off()
normalized_matrix=normalized_matrix[rowMeans(normalized_matrix)>0,]
M <- cor(as.matrix(normalized_matrix),method="spearman")
pdf(file=file.path(output_dir, paste0(matrix_name_flag,"_tumor_normal_UQ_log2(x+1)_BCM_corplot.pdf")),height=7,width=7)
corrplot(M, method = "color",order = "hclust",tl.cex = 0.3,cl.lim=c(min(M),max(M)), is.corr=FALSE)
dev.off()

pdf(file=file.path(output_dir, paste0(matrix_name_flag,"_tumor_normal_raw_log2(x+1)_BCM_pcaplot.pdf")),width=6,height=4)
sample_labels=rep("Tumor",length(colnames(raw_matrix)))
sample_labels[grepl("_A",colnames(raw_matrix))]="Normal"
res.pca <- prcomp(na.omit(t(raw_matrix)),scale = TRUE)
p <- fviz_pca_ind(res.pca, label="none", habillage=sample_labels,addEllipses=TRUE, ellipse.level=0.95)
print(p)
dev.off()

pdf(file=file.path(output_dir, paste0(matrix_name_flag,"_tumor_normal_UQ_log2(x+1)_BCM_pcaplot.pdf")),height=4,width=6)
sample_labels=rep("Tumor",length(colnames(normalized_matrix)))
sample_labels[grepl("_A",colnames(normalized_matrix))]="Normal"
res.pca <- prcomp(na.omit(t(normalized_matrix)),scale = TRUE)
p <- fviz_pca_ind(res.pca, label="none", habillage=sample_labels,addEllipses=TRUE, ellipse.level=0.95)
print(p)
dev.off()

##################################################
##################################################

rm(list=ls())
output_dir <- "RSEM_results_summarized"
matrix_name_flag="gene_rsem_removed_circRNA"
sample_list_for_summary=list.dirs(path = ".", full.names = FALSE, recursive = FALSE)
sample_list_for_summary=sample_list_for_summary[sample_list_for_summary!=output_dir]
example_file_path=file.path(sample_list_for_summary[1], paste0(matrix_name_flag,"_BCM.txt"))
example_table=read.table(example_file_path,header=TRUE,row.names=1,check.name=FALSE,sep="\t")
total_item_number=nrow(example_table)
raw_matrix=matrix(data=0,nrow=total_item_number,ncol=length(sample_list_for_summary))
normalized_matrix=matrix(data=0,nrow=total_item_number,ncol=length(sample_list_for_summary))
rownames(raw_matrix)=rownames(example_table)
colnames(raw_matrix)=sample_list_for_summary

for(i in 1:length(sample_list_for_summary))
{
   input_file_path=file.path(sample_list_for_summary[i], paste0(matrix_name_flag,"_BCM.txt"))
   print(input_file_path)
   sample_table=read.table(input_file_path,header=TRUE,row.names=1,check.name=FALSE,sep="\t")
   raw_matrix[,i]=sample_table[,1]
}
raw_matrix[raw_matrix<0]=0
#normalized_matrix <- normalize.quantiles(as.matrix(raw_matrix), copy = TRUE)
#row.names(normalized_matrix) <- row.names(raw_matrix)
#colnames(normalized_matrix) <- colnames(raw_matrix)

group <- as.factor(colnames(raw_matrix))
cds <- DGEList(counts = raw_matrix, group=group)
cds <- calcNormFactors(cds, method="upperquartile")
scale <- cds$samples$lib.size*cds$samples$norm.factors
normalized_matrix <- round(t(t(raw_matrix)/scale)*mean(scale),2)

write.table(round(log2(raw_matrix+1),2),file=file.path(output_dir, paste0(matrix_name_flag,"_tumor_normal_raw_log2(x+1)_BCM.txt")),quote = FALSE,sep="\t")
write.table(round(log2(normalized_matrix+1),2),file=file.path(output_dir, paste0(matrix_name_flag,"_tumor_normal_UQ_log2(x+1)_BCM.txt")),quote = FALSE,sep="\t")

pdf(file=file.path(output_dir, paste0(matrix_name_flag,"_tumor_normal_raw_log2(x+1)_BCM_boxplot.pdf")),height=7,width=12)
boxplot(log2(raw_matrix+1),las = 2,ylab="Log2(RSEM+1)",main="Distribution of gene expression")
dev.off()
pdf(file=file.path(output_dir, paste0(matrix_name_flag,"_tumor_normal_UQ_log2(x+1)_BCM_boxplot.pdf")),height=7,width=12)
boxplot(log2(normalized_matrix+1),las = 2,ylab="Log2(RSEM+1)",main="Distribution of gene expression")
dev.off()

raw_matrix=raw_matrix[rowMeans(raw_matrix)>0,]
M <- cor(as.matrix(raw_matrix),method="spearman")
pdf(file=file.path(output_dir, paste0(matrix_name_flag,"_tumor_normal_raw_log2(x+1)_BCM_corplot.pdf")),height=7,width=7)
corrplot(M, method = "color",order = "hclust",tl.cex = 0.3,cl.lim=c(min(M),max(M)), is.corr=FALSE)
dev.off()
normalized_matrix=normalized_matrix[rowMeans(normalized_matrix)>0,]
M <- cor(as.matrix(normalized_matrix),method="spearman")
pdf(file=file.path(output_dir, paste0(matrix_name_flag,"_tumor_normal_UQ_log2(x+1)_BCM_corplot.pdf")),height=7,width=7)
corrplot(M, method = "color",order = "hclust",tl.cex = 0.3,cl.lim=c(min(M),max(M)), is.corr=FALSE)
dev.off()

pdf(file=file.path(output_dir, paste0(matrix_name_flag,"_tumor_normal_raw_log2(x+1)_BCM_pcaplot.pdf")),width=6,height=4)
sample_labels=rep("Tumor",length(colnames(raw_matrix)))
sample_labels[grepl("_A",colnames(raw_matrix))]="Normal"
res.pca <- prcomp(na.omit(t(raw_matrix)),scale = TRUE)
p <- fviz_pca_ind(res.pca, label="none", habillage=sample_labels,addEllipses=TRUE, ellipse.level=0.95)
print(p)
dev.off()

pdf(file=file.path(output_dir, paste0(matrix_name_flag,"_tumor_normal_UQ_log2(x+1)_BCM_pcaplot.pdf")),width=6,height=4)
sample_labels=rep("Tumor",length(colnames(raw_matrix)))
sample_labels[grepl("_A",colnames(raw_matrix))]="Normal"
res.pca <- prcomp(na.omit(t(raw_matrix)),scale = TRUE)
p <- fviz_pca_ind(res.pca, label="none", habillage=sample_labels,addEllipses=TRUE, ellipse.level=0.95)
print(p)
dev.off()

##################################################
##################################################
rm(list=ls())
output_dir <- "RSEM_results_summarized"
matrix_name_flag="gene_tpm_with_circRNA"
sample_list_for_summary=list.dirs(path = ".", full.names = FALSE, recursive = FALSE)
sample_list_for_summary=sample_list_for_summary[sample_list_for_summary!=output_dir]
example_file_path=file.path(sample_list_for_summary[1], paste0(matrix_name_flag,"_BCM.txt"))
example_table=read.table(example_file_path,header=TRUE,row.names=1,check.name=FALSE,sep="\t")
total_item_number=nrow(example_table)
raw_matrix=matrix(data=0,nrow=total_item_number,ncol=length(sample_list_for_summary))
normalized_matrix=matrix(data=0,nrow=total_item_number,ncol=length(sample_list_for_summary))
rownames(raw_matrix)=rownames(example_table)
colnames(raw_matrix)=sample_list_for_summary

for(i in 1:length(sample_list_for_summary))
{
   input_file_path=file.path(sample_list_for_summary[i], paste0(matrix_name_flag,"_BCM.txt"))
   print(input_file_path)
   sample_table=read.table(input_file_path,header=TRUE,row.names=1,check.name=FALSE,sep="\t")
   raw_matrix[,i]=sample_table[,1]
}

write.table(round(log2(raw_matrix+1),2),file=file.path(output_dir, paste0(matrix_name_flag,"_tumor_normal_raw_log2(x+1)_BCM.txt")),quote = FALSE,sep="\t")
raw_matrix=raw_matrix[rowMeans(raw_matrix)>0,]
M <- cor(as.matrix(raw_matrix),method="spearman")
pdf(file=file.path(output_dir, paste0(matrix_name_flag,"_tumor_normal_raw_log2(x+1)_BCM_corplot.pdf")),height=7,width=7)
corrplot(M, method = "color",order = "hclust",tl.cex = 0.3,cl.lim=c(min(M),max(M)), is.corr=FALSE)
dev.off()

pdf(file=file.path(output_dir, paste0(matrix_name_flag,"_tumor_normal_raw_log2(x+1)_BCM_pcaplot.pdf")),width=6,height=4)
sample_labels=rep("Tumor",length(colnames(raw_matrix)))
sample_labels[grepl("_A",colnames(raw_matrix))]="Normal"
res.pca <- prcomp(na.omit(t(raw_matrix)),scale = TRUE)
p <- fviz_pca_ind(res.pca, label="none", habillage=sample_labels,addEllipses=TRUE, ellipse.level=0.95)
print(p)
dev.off()

##################################################
##################################################

rm(list=ls())
output_dir <- "RSEM_results_summarized"
matrix_name_flag="gene_tpm_removed_circRNA"
sample_list_for_summary=list.dirs(path = ".", full.names = FALSE, recursive = FALSE)
sample_list_for_summary=sample_list_for_summary[sample_list_for_summary!=output_dir]
example_file_path=file.path(sample_list_for_summary[1], paste0(matrix_name_flag,"_BCM.txt"))
example_table=read.table(example_file_path,header=TRUE,row.names=1,check.name=FALSE,sep="\t")
total_item_number=nrow(example_table)
raw_matrix=matrix(data=0,nrow=total_item_number,ncol=length(sample_list_for_summary))
normalized_matrix=matrix(data=0,nrow=total_item_number,ncol=length(sample_list_for_summary))
rownames(raw_matrix)=rownames(example_table)
colnames(raw_matrix)=sample_list_for_summary

for(i in 1:length(sample_list_for_summary))
{
   input_file_path=file.path(sample_list_for_summary[i], paste0(matrix_name_flag,"_BCM.txt"))
   print(input_file_path)
   sample_table=read.table(input_file_path,header=TRUE,row.names=1,check.name=FALSE,sep="\t")
   raw_matrix[,i]=sample_table[,1]
}

write.table(round(log2(raw_matrix+1),2),file=file.path(output_dir, paste0(matrix_name_flag,"_tumor_normal_raw_log2(x+1)_BCM.txt")),quote = FALSE,sep="\t")

raw_matrix=raw_matrix[rowMeans(raw_matrix)>0,]
M <- cor(as.matrix(raw_matrix),method="spearman")
pdf(file=file.path(output_dir, paste0(matrix_name_flag,"_tumor_normal_raw_log2(x+1)_BCM_corplot.pdf")),height=7,width=7)
corrplot(M, method = "color",order = "hclust",tl.cex = 0.3,cl.lim=c(min(M),max(M)), is.corr=FALSE)
dev.off()

pdf(file=file.path(output_dir, paste0(matrix_name_flag,"_tumor_normal_raw_log2(x+1)_BCM_pcaplot.pdf")),width=6,height=4)
sample_labels=rep("Tumor",length(colnames(raw_matrix)))
sample_labels[grepl("_A",colnames(raw_matrix))]="Normal"
res.pca <- prcomp(na.omit(t(raw_matrix)),scale = TRUE)
p <- fviz_pca_ind(res.pca, label="none", habillage=sample_labels,addEllipses=TRUE, ellipse.level=0.95)
print(p)
dev.off()

##################################################
##################################################

rm(list=ls())
output_dir <- "RSEM_results_summarized"
matrix_name_flag="gene_fpkm_with_circRNA"
sample_list_for_summary=list.dirs(path = ".", full.names = FALSE, recursive = FALSE)
sample_list_for_summary=sample_list_for_summary[sample_list_for_summary!=output_dir]
example_file_path=file.path(sample_list_for_summary[1], paste0(matrix_name_flag,"_BCM.txt"))
example_table=read.table(example_file_path,header=TRUE,row.names=1,check.name=FALSE,sep="\t")
total_item_number=nrow(example_table)
raw_matrix=matrix(data=0,nrow=total_item_number,ncol=length(sample_list_for_summary))
normalized_matrix=matrix(data=0,nrow=total_item_number,ncol=length(sample_list_for_summary))
rownames(raw_matrix)=rownames(example_table)
colnames(raw_matrix)=sample_list_for_summary

for(i in 1:length(sample_list_for_summary))
{
   input_file_path=file.path(sample_list_for_summary[i], paste0(matrix_name_flag,"_BCM.txt"))
   print(input_file_path)
   sample_table=read.table(input_file_path,header=TRUE,row.names=1,check.name=FALSE,sep="\t")
   raw_matrix[,i]=sample_table[,1]
}

write.table(round(log2(raw_matrix+1),2),file=file.path(output_dir, paste0(matrix_name_flag,"_tumor_normal_raw_log2(x+1)_BCM.txt")),quote = FALSE,sep="\t")

raw_matrix=raw_matrix[rowMeans(raw_matrix)>0,]
M <- cor(as.matrix(raw_matrix),method="spearman")
pdf(file=file.path(output_dir, paste0(matrix_name_flag,"_tumor_normal_raw_log2(x+1)_BCM_corplot.pdf")),height=7,width=7)
corrplot(M, method = "color",order = "hclust",tl.cex = 0.3,cl.lim=c(min(M),max(M)), is.corr=FALSE)
dev.off()

pdf(file=file.path(output_dir, paste0(matrix_name_flag,"_tumor_normal_raw_log2(x+1)_BCM_pcaplot.pdf")),width=6,height=4)
sample_labels=rep("Tumor",length(colnames(raw_matrix)))
sample_labels[grepl("_A",colnames(raw_matrix))]="Normal"
res.pca <- prcomp(na.omit(t(raw_matrix)),scale = TRUE)
p <- fviz_pca_ind(res.pca, label="none", habillage=sample_labels,addEllipses=TRUE, ellipse.level=0.95)
print(p)
dev.off()

##################################################
##################################################

rm(list=ls())
output_dir <- "RSEM_results_summarized"
matrix_name_flag="gene_fpkm_removed_circRNA"
sample_list_for_summary=list.dirs(path = ".", full.names = FALSE, recursive = FALSE)
sample_list_for_summary=sample_list_for_summary[sample_list_for_summary!=output_dir]
example_file_path=file.path(sample_list_for_summary[1], paste0(matrix_name_flag,"_BCM.txt"))
example_table=read.table(example_file_path,header=TRUE,row.names=1,check.name=FALSE,sep="\t")
total_item_number=nrow(example_table)
raw_matrix=matrix(data=0,nrow=total_item_number,ncol=length(sample_list_for_summary))
normalized_matrix=matrix(data=0,nrow=total_item_number,ncol=length(sample_list_for_summary))
rownames(raw_matrix)=rownames(example_table)
colnames(raw_matrix)=sample_list_for_summary

for(i in 1:length(sample_list_for_summary))
{
   input_file_path=file.path(sample_list_for_summary[i], paste0(matrix_name_flag,"_BCM.txt"))
   print(input_file_path)
   sample_table=read.table(input_file_path,header=TRUE,row.names=1,check.name=FALSE,sep="\t")
   raw_matrix[,i]=sample_table[,1]
}

write.table(round(log2(raw_matrix+1),2),file=file.path(output_dir, paste0(matrix_name_flag,"_tumor_normal_raw_log2(x+1)_BCM.txt")),quote = FALSE,sep="\t")

raw_matrix=raw_matrix[rowMeans(raw_matrix)>0,]
M <- cor(as.matrix(raw_matrix),method="spearman")
pdf(file=file.path(output_dir, paste0(matrix_name_flag,"_tumor_normal_raw_log2(x+1)_BCM_corplot.pdf")),height=7,width=7)
corrplot(M, method = "color",order = "hclust",tl.cex = 0.3,cl.lim=c(min(M),max(M)), is.corr=FALSE)
dev.off()

pdf(file=file.path(output_dir, paste0(matrix_name_flag,"_tumor_normal_raw_log2(x+1)_BCM_pcaplot.pdf")),width=6,height=4)
sample_labels=rep("Tumor",length(colnames(raw_matrix)))
sample_labels[grepl("_A",colnames(raw_matrix))]="Normal"
res.pca <- prcomp(na.omit(t(raw_matrix)),scale = TRUE)
p <- fviz_pca_ind(res.pca, label="none", habillage=sample_labels,addEllipses=TRUE, ellipse.level=0.95)
print(p)
dev.off()

##################################################
##################################################

rm(list=ls())
output_dir <- "RSEM_results_summarized"
matrix_name_flag="linear_isoform_rsem"
sample_list_for_summary=list.dirs(path = ".", full.names = FALSE, recursive = FALSE)
sample_list_for_summary=sample_list_for_summary[sample_list_for_summary!=output_dir]
example_file_path=file.path(sample_list_for_summary[1], paste0(matrix_name_flag,"_BCM.txt"))
example_table=read.table(example_file_path,header=TRUE,row.names=1,check.name=FALSE,sep="\t")
total_item_number=nrow(example_table)
raw_matrix=matrix(data=0,nrow=total_item_number,ncol=length(sample_list_for_summary))
normalized_matrix=matrix(data=0,nrow=total_item_number,ncol=length(sample_list_for_summary))
rownames(raw_matrix)=rownames(example_table)
colnames(raw_matrix)=sample_list_for_summary

for(i in 1:length(sample_list_for_summary))
{
   input_file_path=file.path(sample_list_for_summary[i], paste0(matrix_name_flag,"_BCM.txt"))
   print(input_file_path)
   sample_table=read.table(input_file_path,header=TRUE,row.names=1,check.name=FALSE,sep="\t")
   raw_matrix[,i]=sample_table[,1]
}

write.table(round(log2(raw_matrix+1),2),file=file.path(output_dir, paste0(matrix_name_flag,"_tumor_normal_raw_log2(x+1)_BCM.txt")),quote = FALSE,sep="\t")

##################################################
##################################################

rm(list=ls())
output_dir <- "RSEM_results_summarized"
matrix_name_flag="linear_isoform_tpm"
sample_list_for_summary=list.dirs(path = ".", full.names = FALSE, recursive = FALSE)
sample_list_for_summary=sample_list_for_summary[sample_list_for_summary!=output_dir]
example_file_path=file.path(sample_list_for_summary[1], paste0(matrix_name_flag,"_BCM.txt"))
example_table=read.table(example_file_path,header=TRUE,row.names=1,check.name=FALSE,sep="\t")
total_item_number=nrow(example_table)
raw_matrix=matrix(data=0,nrow=total_item_number,ncol=length(sample_list_for_summary))
normalized_matrix=matrix(data=0,nrow=total_item_number,ncol=length(sample_list_for_summary))
rownames(raw_matrix)=rownames(example_table)
colnames(raw_matrix)=sample_list_for_summary

for(i in 1:length(sample_list_for_summary))
{
   input_file_path=file.path(sample_list_for_summary[i], paste0(matrix_name_flag,"_BCM.txt"))
   print(input_file_path)
   sample_table=read.table(input_file_path,header=TRUE,row.names=1,check.name=FALSE,sep="\t")
   raw_matrix[,i]=sample_table[,1]
}

write.table(round(log2(raw_matrix+1),2),file=file.path(output_dir, paste0(matrix_name_flag,"_tumor_normal_raw_log2(x+1)_BCM.txt")),quote = FALSE,sep="\t")

##################################################
##################################################

rm(list=ls())
output_dir <- "RSEM_results_summarized"
matrix_name_flag="linear_isoform_fpkm"
sample_list_for_summary=list.dirs(path = ".", full.names = FALSE, recursive = FALSE)
sample_list_for_summary=sample_list_for_summary[sample_list_for_summary!=output_dir]
example_file_path=file.path(sample_list_for_summary[1], paste0(matrix_name_flag,"_BCM.txt"))
example_table=read.table(example_file_path,header=TRUE,row.names=1,check.name=FALSE,sep="\t")
total_item_number=nrow(example_table)
raw_matrix=matrix(data=0,nrow=total_item_number,ncol=length(sample_list_for_summary))
normalized_matrix=matrix(data=0,nrow=total_item_number,ncol=length(sample_list_for_summary))
rownames(raw_matrix)=rownames(example_table)
colnames(raw_matrix)=sample_list_for_summary

for(i in 1:length(sample_list_for_summary))
{
   input_file_path=file.path(sample_list_for_summary[i], paste0(matrix_name_flag,"_BCM.txt"))
   print(input_file_path)
   sample_table=read.table(input_file_path,header=TRUE,row.names=1,check.name=FALSE,sep="\t")
   raw_matrix[,i]=sample_table[,1]
}

write.table(round(log2(raw_matrix+1),2),file=file.path(output_dir, paste0(matrix_name_flag,"_tumor_normal_raw_log2(x+1)_BCM.txt")),quote = FALSE,sep="\t")

##################################################
##################################################

rm(list=ls())
output_dir <- "RSEM_results_summarized"
matrix_name_flag="circRNA_rsem"
sample_list_for_summary=list.dirs(path = ".", full.names = FALSE, recursive = FALSE)
sample_list_for_summary=sample_list_for_summary[sample_list_for_summary!=output_dir]
circ_RNA_list=rep("NA",1)
for(i in 1:length(sample_list_for_summary))
{
   input_file_path=file.path(sample_list_for_summary[i], paste0(matrix_name_flag,"_BCM.txt"))
   print(input_file_path)
   sample_table=read.table(input_file_path,header=TRUE,row.names=1,check.name=FALSE,sep="\t")
   circ_RNA_list_sample=rownames(sample_table)
   circ_RNA_list=c(circ_RNA_list,circ_RNA_list_sample)
}
circ_RNA_list=unique(circ_RNA_list)
circ_RNA_list=circ_RNA_list[-1]

total_item_number=length(circ_RNA_list)
raw_matrix=matrix(data=0,nrow=total_item_number,ncol=length(sample_list_for_summary))
rownames(raw_matrix)=circ_RNA_list
colnames(raw_matrix)=sample_list_for_summary

for(i in 1:length(sample_list_for_summary))
{
   input_file_path=file.path(sample_list_for_summary[i], paste0(matrix_name_flag,"_BCM.txt"))
   print(input_file_path)
   sample_table=read.table(input_file_path,header=TRUE,row.names=1,check.name=FALSE,sep="\t")
   sample_table[setdiff(circ_RNA_list,rownames(sample_table)),]=0
   sample_table[,2]=sample_table[,1]
   sample_table=as.data.frame(sample_table[circ_RNA_list,])
   raw_matrix[,i]=sample_table[,1]
}

write.table(round(log2(raw_matrix+1),2),file=file.path(output_dir, paste0(matrix_name_flag,"_tumor_normal_raw_log2(x+1)_BCM.txt")),quote = FALSE,sep="\t")

##################################################
##################################################

rm(list=ls())
output_dir <- "RSEM_results_summarized"
matrix_name_flag="circRNA_fpkm"
sample_list_for_summary=list.dirs(path = ".", full.names = FALSE, recursive = FALSE)
sample_list_for_summary=sample_list_for_summary[sample_list_for_summary!=output_dir]
circ_RNA_list=rep("NA",1)
for(i in 1:length(sample_list_for_summary))
{
   input_file_path=file.path(sample_list_for_summary[i], paste0(matrix_name_flag,"_BCM.txt"))
   print(input_file_path)
   sample_table=read.table(input_file_path,header=TRUE,row.names=1,check.name=FALSE,sep="\t")
   circ_RNA_list_sample=rownames(sample_table)
   circ_RNA_list=c(circ_RNA_list,circ_RNA_list_sample)
}
circ_RNA_list=unique(circ_RNA_list)
circ_RNA_list=circ_RNA_list[-1]

total_item_number=length(circ_RNA_list)
raw_matrix=matrix(data=0,nrow=total_item_number,ncol=length(sample_list_for_summary))
rownames(raw_matrix)=circ_RNA_list
colnames(raw_matrix)=sample_list_for_summary

for(i in 1:length(sample_list_for_summary))
{
   input_file_path=file.path(sample_list_for_summary[i], paste0(matrix_name_flag,"_BCM.txt"))
   print(input_file_path)
   sample_table=read.table(input_file_path,header=TRUE,row.names=1,check.name=FALSE,sep="\t")
   sample_table[setdiff(circ_RNA_list,rownames(sample_table)),]=0
   sample_table[,2]=sample_table[,1]
   sample_table=as.data.frame(sample_table[circ_RNA_list,])
   raw_matrix[,i]=sample_table[,1]
}

write.table(round(log2(raw_matrix+1),2),file=file.path(output_dir, paste0(matrix_name_flag,"_tumor_normal_raw_log2(x+1)_BCM.txt")),quote = FALSE,sep="\t")

##################################################
##################################################

rm(list=ls())
output_dir <- "RSEM_results_summarized"
matrix_name_flag="circRNA_tpm"
sample_list_for_summary=list.dirs(path = ".", full.names = FALSE, recursive = FALSE)
sample_list_for_summary=sample_list_for_summary[sample_list_for_summary!=output_dir]

circ_RNA_list=rep("NA",1)
for(i in 1:length(sample_list_for_summary))
{
   input_file_path=file.path(sample_list_for_summary[i], paste0(matrix_name_flag,"_BCM.txt"))
   print(input_file_path)
   sample_table=read.table(input_file_path,header=TRUE,row.names=1,check.name=FALSE,sep="\t")
   circ_RNA_list_sample=rownames(sample_table)
   circ_RNA_list=c(circ_RNA_list,circ_RNA_list_sample)
}
circ_RNA_list=unique(circ_RNA_list)
circ_RNA_list=circ_RNA_list[-1]

total_item_number=length(circ_RNA_list)
raw_matrix=matrix(data=0,nrow=total_item_number,ncol=length(sample_list_for_summary))
rownames(raw_matrix)=circ_RNA_list
colnames(raw_matrix)=sample_list_for_summary

for(i in 1:length(sample_list_for_summary))
{
   input_file_path=file.path(sample_list_for_summary[i], paste0(matrix_name_flag,"_BCM.txt"))
   print(input_file_path)
   sample_table=read.table(input_file_path,header=TRUE,row.names=1,check.name=FALSE,sep="\t")
   sample_table[setdiff(circ_RNA_list,rownames(sample_table)),]=0
   sample_table[,2]=sample_table[,1]
   sample_table=as.data.frame(sample_table[circ_RNA_list,])
   raw_matrix[,i]=sample_table[,1]
}

write.table(round(log2(raw_matrix+1),2),file=file.path(output_dir, paste0(matrix_name_flag,"_tumor_normal_raw_log2(x+1)_BCM.txt")),quote = FALSE,sep="\t")


##################################################
##################################################

rm(list=ls())
output_dir <- "RSEM_results_summarized"
matrix_name_flag="circRNA_rsem"
sample_list_for_summary=list.dirs(path = ".", full.names = FALSE, recursive = FALSE)
sample_list_for_summary=sample_list_for_summary[sample_list_for_summary!=output_dir]
circ_RNA_list=rep("NA",1)
for(i in 1:length(sample_list_for_summary))
{
   input_file_path=file.path(sample_list_for_summary[i], paste0(matrix_name_flag,"_BCM.txt"))
   print(input_file_path)
   sample_table=read.table(input_file_path,header=TRUE,row.names=1,check.name=FALSE,sep="\t")
   circ_RNA_list_sample=rownames(sample_table)
   circ_RNA_list=c(circ_RNA_list,circ_RNA_list_sample)
}
circ_RNA_list=unique(circ_RNA_list)
circ_RNA_list=circ_RNA_list[-1]

total_item_number=length(circ_RNA_list)
raw_matrix=matrix(data=0,nrow=total_item_number,ncol=length(sample_list_for_summary))
rownames(raw_matrix)=circ_RNA_list
colnames(raw_matrix)=sample_list_for_summary

for(i in 1:length(sample_list_for_summary))
{
   input_file_path=file.path(sample_list_for_summary[i], paste0(matrix_name_flag,"_BCM.txt"))
   print(input_file_path)
   sample_table=read.table(input_file_path,header=TRUE,row.names=1,check.name=FALSE,sep="\t")
   sample_table[setdiff(circ_RNA_list,rownames(sample_table)),]=0
   sample_table[,2]=sample_table[,1]
   sample_table=as.data.frame(sample_table[circ_RNA_list,])
   raw_matrix[,i]=sample_table[,1]
}

linear_gene_expression=read.table(file.path(output_dir, "gene_fpkm_removed_circRNA_tumor_normal_raw_log2(x+1)_BCM.txt"),header=TRUE,row.names=1,check.name=FALSE,sep="\t")

combined_expression=rbind(linear_gene_expression,raw_matrix)
combined_expression[combined_expression<0]=0
group <- as.factor(colnames(combined_expression))
cds <- DGEList(counts = combined_expression, group=group)
cds <- calcNormFactors(cds, method="upperquartile")
scale <- cds$samples$lib.size*cds$samples$norm.factors
normalized_matrix <- round(t(t(combined_expression)/scale)*mean(scale),2)
normalized_matrix=normalized_matrix[grepl("circ_",rownames(normalized_matrix)),]
write.table(round(log2(normalized_matrix+1),2),file=file.path(output_dir, paste0(matrix_name_flag,"_tumor_normal_UQ_log2(x+1)_BCM.txt")),quote = FALSE,sep="\t")
