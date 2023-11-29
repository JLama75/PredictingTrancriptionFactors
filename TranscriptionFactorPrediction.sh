#Model 2

#ATAC-seq
cat Pachy.narrowPeak | awk 'BEGIN{OFS="\t"};{print $1,$2,$3}' | sort -k1,1 -k2,2n > Pachy.narrowPeak.sorted.bed #58,394
wc -l Pachy.narrowPeak
wc -l Pachy.narrowPeak.sorted.bed

#dREG 
zcat P_merged.dREG.peak.score.bed.gz | awk 'BEGIN{OFS="\t"};{print $1,$2,$3}' | sort -k1,1 -k2,2n > P_dREG_sorted.bed #65,533

#For experimental: (Pachy open chromatin - Pachy dREGs) in sex chromosomes
#Only report those entries in A that have no overlap in B.
bedtools intersect -a Pachy.narrowPeak.sorted.bed -b P_dREG_sorted.bed -v > P_open_withoutdREG.bed #41527
cat P_open_withoutdREG.bed | awk '{print$1}' | uniq -c > P_open_withoutdREG.counts

#subset only in X and Y chromosome #150 +8
cat P_open_withoutdREG.bed | awk 'BEGIN{OFS="\t"} ($1=="chrX" || $1=="chrY") {print $0}' > chrXY_P_open_withoutdREG.bed #150 + 8
cat P_open_withoutdREG.bed | awk '{print$1}' | uniq -c > chrXY_open_withoutdREG.counts

#For Background set 1 (open chromatin with dREG sites + dREG sites)
bedtools intersect -a Pachy.narrowPeak.sorted.bed -b P_dREG_sorted.bed > P_open_withdREG.bed #34,182
#merge the open_withdREG with Pachy_dREGs
cat P_open_withdREG.bed P_dREG_sorted.bed | sort -k1,1 -k2,2n > tmp.bed #99715
bedtools merge -i tmp.bed > P_open_withdREG_dREG.bed #65,533
rm tmp.bed

#For Background set 2
#open chromatin with dREGs + all pachy dREGs + open chromatin without dREGs in autosomes
cat P_open_withoutdREG.bed | grep -v "chrX" | grep -v "chrY" > Autosomes_P_open_withoutdREG.bed #41527-150-8=41369
cat Autosomes_P_open_withoutdREG.bed P_open_withdREG_dREG.bed | sort -k1,1 -k2,2n > tmp.bed #106910
bedtools merge -i tmp.bed > Autosomes_open_dREG_merge.bed #106,905

#Generate uniq ids for HOMER
awk '{printf "peak_%d\t%s\n", NR, $0}' Autosomes_open_dREG_merge.bed > background2.bed
awk '{printf "peak_%d\t%s\n", NR, $0}' P_open_withdREG_dREG.bed > P_open_withdREG_dREG2.bed
awk '{printf "peak_%d\t%s\n", NR, $0}' chrXY_P_open_withoutdREG.bed > chrXY_P_open_withoutdREG2.bed

cat background2.bed | awk 'BEGIN{OFS="\t"}; {print $2,$3,$4,$1}' > background_2.bed
awk 'BEGIN{OFS="\t"}; {print $2,$3,$4,$1}' P_open_withdREG_dREG2.bed > P_open_withdREG_dREG.bed
awk 'BEGIN{OFS="\t"}; {print $2,$3,$4,$1}' chrXY_P_open_withoutdREG2.bed > chrXY_P_open_withoutdREG.bed

rm background2.bed P_open_withdREG_dREG2.bed chrXY_P_open_withoutdREG2.bed Autosomes_open_dREG_merge.bed


#install homer in /workdir/jl3285/HOMER/homer
#mkdir homer
#cd homer
#wget http://homer.ucsd.edu/homer/configureHomer.pl
#perl configureHomer.pl -install homer

#install the database of your species
#perl configureHomer.pl -install mm10

source $HOME/miniconda3/bin/activate renv
export PATH=/workdir/jl3285/HOMER/homer/bin:$PATH

findMotifsGenome.pl chrXY_P_open_withoutdREG.bed mm10 P_model_2a_b1_MotifOutput -bg P_open_withdREG_dREG.bed >& bed.log &
findMotifsGenome.pl chrXY_P_open_withoutdREG.bed mm10 P_model_2a_b2_MotifOutput -bg background_2.bed  >& motif.log &
findMotifsGenome.pl chrXY_P_open_withoutdREG.bed mm10 P_model_2a_b2_MotifOutput >& motif.log & #random generated background
