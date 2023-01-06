#Pachytene spermatocyte vs. Leptotene/Zygotene spermatocyte

library(dplyr)
df1 <- read.table("P.vs.LZ.gene.deseq.txt", header=TRUE)
df <- df1[,c(1,2,3,4,5,6,7,8,12)]
sum(df$gene.padj < 0.05, na.rm=TRUE) #14569
df.padj <- df %>% filter(gene.padj < 0.05) %>% arrange(desc(gene.padj))

df.padj.l2fc.up <- df.padj %>% arrange(desc(gene.log2FoldChange)) %>% filter(gene.log2FoldChange > 0.5) %>% arrange(desc(gene.log2FoldChange)) #upregulated genes in Pachytene
sum(df.padj$gene.log2FoldChange > 0.5, na.rm=TRUE) #7439
sum(df.padj$gene.log2FoldChange > 0, na.rm=TRUE) #7601


df.padj.l2fc.down <- df.padj %>% filter(gene.log2FoldChange < -0.5) %>% arrange(desc(gene.log2FoldChange)) #upregulated genes in Pachytene #6793
sum(df.padj$gene.log2FoldChange < -0.5, na.rm=TRUE)#6793
sum(df.padj$gene.log2FoldChange < 0, na.rm=TRUE) #6968

test <- df %>% filter(gene.log2FoldChange > 0) %>% arrange(desc(gene.log2FoldChange))
head(test) #look at the max l2fc value. does it match the df.padj.l2fc.up

write.table(df.padj.l2fc.up, "P.vs.LZ.DEG_padj_l2fc_upregulated.txt",sep="\t")
write.table(df.padj.l2fc.down, "P.vs.LZ.DEG_padj_l2fc_downregulated.txt",sep="\t")
