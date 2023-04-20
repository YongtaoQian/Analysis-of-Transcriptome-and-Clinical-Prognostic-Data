library("EnhancedVolcano")
args = commandArgs(trailingOnly=TRUE)
allDiff = read.table(args[1], sep = "\t", header = TRUE)

res <- allDiff
pdf(file="Volcano.pdf",onefile = FALSE,width = 12.5,height =9)
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'logFC',
                y = 'PValue',
                title = 'Retroperitoneal Liposarcoma',
                pCutoff = 0.05,
                FCcutoff = 1.0,
                pointSize = 1.0,
                labSize = 6.0,
                legendPosition = 'top',
                legendLabSize = 16,
                legendIconSize = 5.0,
                drawConnectors = F,
                gridlines.major = FALSE,     # 背景网格
                gridlines.minor = FALSE,
                colAlpha = 1)
dev.off()

library(readxl)
library(clusterProfiler)
library(ggplot2)
library(enrichplot)
library(GOplot)
library(DOSE)
library(org.Hs.eg.db)

data <- read.table("diff.xls",sep="\t")  ##另存一下即可
eg <- bitr(rownames(data), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
gene <- eg$ENTREZID
epath <- enrichKEGG(gene = gene,
                    organism ="hsa",
                    keyType = "kegg",
                    pvalueCutoff = 1,
                    qvalueCutoff = 1,
                    use_internal_data = F)
epath <- setReadable(epath,OrgDb = org.Hs.eg.db,keyType = "ENTREZID" ) 
write.csv(as.data.frame(epath@result), file="Path_enrichment1.csv")
p1=barplot(epath, showCategory=19,title="Path_enrichment",label_format=100)#条状图，按p从小到大排的
p2=dotplot(epath,showCategory=19,title="Path_enrichment",label_format=100)#点图，按富集的数从大到小的
ggsave(p1,filename ="Path_enrichment_bar1.pdf",width =10,height = 6 )
ggsave(p2,filename ="Path_enrichment_dot1.pdf",width =10,height = 6 )

