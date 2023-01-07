
#Find overlap regions of dREG sites (Predicted Transcriptional Regulatory Elements) and downregulated genes

cat P.vs.LZ.DEG_padj_l2fc_downregulated.txt | sed -e s/\"//g > downregulated_gene_pos.bed 
awk '{if (NR!=1){print $2, $3, $4, $7}}' downregulated_gene_pos.bed | sort -k1,1 -k2,2n > downregulated_gene_pos_sorted.bed #6793

#ready the dREG files
awk '{print $1, $2, $3}' P_dREG.peak.score.bed | sort -k1,1 -k2,2n > P_dREG.bed #65533
#1000 bp near genes
awk '{print $1, $2, $3+1000}' P_dREG.peak.score.bed | sort -k1,1 -k2,2n > P_dREG.neargene.bed #65533


#ensure the files are tab delimited for bedtools
source $HOME/miniconda3/bin/activate/renv
#Using R/2.2.0 version 

R
df1 <- read.table("P_dREG.neargene.bed", header=F)
df2 <- read.table("P_dREG.bed", header=F)
df3 <- read.table("downregulated_gene_pos_sorted.bed", header=F)
head(df1)
colnames(df1) <- NULL
colnames(df2) <- NULL
colnames(df3) <- NULL
df1[1,3]
         
write.table(df1,"P_dREG.neargene.bed", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
write.table(df2,"P_dREG.bed", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
write.table(df3,"downregulated_gene_pos_sorted.bed", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

#selecting dREG sites that overlap with downregulated genes
#Use bedtools window: reports all -a that are 1000 bp upstream of -b. define upstream and downstream depending on strand
#report all the genes that are within 1000bp downstream of dREG sites. adds 1000 bp to chrom end of -a (dREG)files.

bedtools window -a P_dREG.bed  -b downregulated_gene_pos_sorted.bed -l 0 -r 1000 -sw > P.vs.LZ.dREG_near_downregulated_2.bed # 15770
bedtools window -a P_dREG.bed  -b downregulated_gene_pos_sorted.bed -l 0 -r 1000  | wc -l
bedtools window -a P_dREG.bed  -b downregulated_gene_pos_sorted.bed -l 0 -r 0  | wc -l
#there are some repeated dREG sites due to it being overlapped with multiple genes or of different strands

#so counting the unique dREGs
awk '{print $1, $2, $3}' P.vs.LZ.dREG_near_downregulated_2.bed | sort -u -k1,1 -k2,2n -k3,3n -s > P.vs.LZ.dREG_near_downregulated_uniq_2.bed
wc -l P.vs.LZ.dREG_near_downregulated_uniq_2.bed #14376 -->(14376/65533*100)--> 21.93% #this is the file we want!

#if you want no repeated dREG sites to be reported use -u instead:
bedtools window -a P_dREG.bed  -b downregulated_gene_pos_sorted2.bed -l 0 -r 1000 -sw -u | wc -l #14376

##for number of hits
bedtools window -a P_dREG.bed  -b downregulated_gene_pos_sorted2.bed -l 0 -r 1000 -sw -c > test2.bed #65533
