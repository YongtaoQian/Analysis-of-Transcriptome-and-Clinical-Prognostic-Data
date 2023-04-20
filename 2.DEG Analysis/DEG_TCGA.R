library(edgeR)
library(limma)
args = commandArgs(trailingOnly=TRUE)
#将第一列换成行名
DATA = read.table(args[1], sep = "\t", header = TRUE,row.names=1)

#分组信息
group_list <- ifelse(as.numeric(substr(colnames(DATA),14,15))<10,"tumor","normal")
group_list <- factor(group_list,levels = c("normal","tumor"))
table(group_list)
dge <- DGEList(counts=DATA,group=group_list)
#差异分析
dge$samples$lib.size <- colSums(dge$counts)
dge <- calcNormFactors(dge)
design <- model.matrix(~0+group_list)
rownames(design)<-colnames(dge)
colnames(design)<-levels(group_list)
dge <- estimateGLMCommonDisp(dge,design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)
fit <- glmFit(dge, design)
fit2 <- glmLRT(fit, contrast=c(-1,1))
DEG2=topTags(fit2, n=nrow(DATA))
DEG2=as.data.frame(DEG2)
write.table(DEG2,file="allgene.xls",sep="\t",quote=F)

diffSig <- DEG2[with(DEG2, (abs(logFC)>1 & PValue < 0.05)), ]
write.table(diffSig,file="diff.xls",sep="\t",quote=F)
diffUp <- DEG2[with(DEG2, (logFC>1 & PValue < 0.05)), ]
write.table(diffUp,file="up.xls",sep="\t",quote=F)
diffDown <- DEG2[with(DEG2, (logFC<(1) & PValue < 0.05)), ]
write.table(diffDown,file="down.xls",sep="\t",quote=F)
