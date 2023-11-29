scratch=/workdir/jl3285/ATAC-seq

final=/local/storage/projects/sexChromosome_dREG_ATAC/workdir_copy/ATAC-seq

echo $scratch

cp *.gz -t $scratch

cd $scratch

rm *_merge*.gz
rm ATAC_D_merge.filtered*
#ATAC_Pachy_05-19-17.filtered.sorted.nd.bam.gz
echo -e "check\n" > out

#starting with the bam files....
for file in *.gz
do

NAMES=$(echo "$file" | cut -d . -f 1 | cut -d _ -f 1-3 | sort | uniq)    
echo processing $file #ATAC_Pachy_05-19-17.filtered.sorted.nd.bam.gz
echo processed the names into "$NAMES" #ATAC_Pachy_05-19-17

echo checking if reads were filtered to remove thoes of low mapping quality...
echo -e "\nCheck if reads for ${NAMES}.filtered.sorted.nd.bam.gz were filtered to remove thoes of low mapping quality..." >> out
zcat {$NAMES}.filtered.sorted.nd.bam.gz | samtools view | awk '($5<20) {print $0}' | wc -l >> out

gzip -d {$NAMES}.filtered.sorted.nd.bam.gz 

#zcat ${NAMES}.filtered.sorted.nd.bam.gz | samtools view -b -q 10 | samtools sort > ${NAMES}.filtered.sorted.qc.nd.bam.gz 
#filtered=mapq>20
#sorted=sorted
#qc=mapq>20
#nd=no PCR duplicates

echo sorting by name..
cat ${NAMES}.filtered.sorted.nd.bam | samtools view -b | samtools sort -n -@ 10 > ${NAMES}.filtered.sorted2.nd.bam

done


#source /home/jl3285/miniconda3/bin/activate Genrich
#echo Genrich activated

for file in `cat text.txt`
do

NAMES=$(echo "$file" | cut -d _ -f 2 | sort | uniq)
echo processing $file #ATAC_Pachy_05-19-17.filtered.sorted.nd.bam.gz
echo processed the names into "$NAMES" #Pachy

/workdir/jl3285/ATAC-seq/Genrich/Genrich -t ATAC_${NAMES}_05-19-17.filtered.sorted2.nd.bam,ATAC_${NAMES}_05-31-17.filtered.sorted2.nd.bam -o ${NAMES}.narrowPeak \
-j -f ${NAMES}.log -q 0.05 -a 20.0 -v -e chrM -y -r

done
 #-r remove duplicates, also accounts for multi-mapping reads which is not accounted by PICARD #optional
 #-j use ATAC-seq mode
 #-f ouput bedgraph-ish file for p/q values (P-value and FDR corrected p-value)
 #-q max q-value (FDR-adjusted, defualt=1)
 #-a max AUC (area under the curve) for a peak (default = 200.0)
 #-v Option to print status updates/counts to stderr
 #-e Comma-separated list of chromosomes to exclude
 #-y Keep unpaired alignments (def. false) 
 
cp *narrowPeak *log *out -t $final
