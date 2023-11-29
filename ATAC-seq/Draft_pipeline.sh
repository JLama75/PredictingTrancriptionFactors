#!/usr/bin/env bash
FASTQDIR="."
# genome file (bowtie2 indexed)
GENOME="/local/storage/data/short_read_index/hg19/bowtie2/hg19"
# Remove the Nextera R2 sequence from ends of R1
ADSEQ_R1="CTGTCTCTTATACACATCT"
ADSEQ_R2="AGATGTGTATAAGAGACAG"
# File containing hg19 chromosome lengths
CHRSZ="/hg19.chromInfo.txt"
# number of processors to use
nproc=32

#First do fastqc 
#for i in ${FASTQDIR}/*.fastq*
#do
#	fpath=${i%.fastq*}
#	[ ! -d ${fpath}_fastqc ] && fastqc -f fastq -t $nproc --casava ${fpath}*.fastq.gz;
#-t --> threads
# -casava --> Files  come  from raw casava output. Files in the same sample group (differing only by the group number) will be analysed as a set rather than individually.
#done

# Process each FASTQ file in the specified directory
for i in ${FASTQDIR}/*_R1.fastq.gz
do
	fpath=${i%_R1.fastq.gz} #${variable%[pattern]} will remove [pattern] from ${variable} if it is at the end
	fname=${fpath##*/} #${variable##[pattern]} will remove [pattern] from ${variable} if it is in the beginning
	echo $fname

	if [ ! -e $fpath.bam ]
	then

		##1. Trim Nextera adapters
		~/.local/bin/cutadapt -a $ADSEQ_R1 -A $ADSEQ_R2 ${fpath}_R1.fastq.gz ${fpath}_R2.fastq.gz -o ${fpath}_trimmed_R1.fastq.gz -p ${fpath}_trimmed_R2.fastq.gz
		## 2. Align reads to hg19 using bowtie2
		# --seed = set random seed (for reproducibility)
		# --no-mixed = don't allow unpaired alignments
		# --no-discordant = remove pairs that are not near each other (default 500 bp)
		# --no-dovetail = remove overlapping alignment pairs (too short)
		bowtie2 -p $nproc -x $GENOME --seed 1026 \
		--no-discordant --no-dovetail --no-unal --no-mixed \
		-1 ${fpath}_trimmed_R1.fastq.gz -2 ${fpath}_trimmed_R2.fastq.gz |
		samtools view -buS -q 10 -F 256 '-' | #-q --> skip alignments with MAPQ smaller than 10.
		 # -F --> exclude Flags with 256 (secondary alignments) maybe for unique alignments
		## 3. Sort & compress in BAM format
		samtools sort -@ $nproc '-' > $fpath.bam;

		# exit with error if BAM wasn't made
		[ -e $fpath.bam ] || exit 1;
	fi
    
    if [ ! -e $fpath.bedGraph ]
    then
        # 4. Compute genome coverage and write bedGraph files
        # -bg = output as bedgraph
        # -fs = use paired-end fragment size, instead of read lengths
        # -ibam = input is BAM file
        # -i = specify input
        # -g = specify chrom lengths
        bedtools genomecov -ibam $fpath.bam -pc -bg > $fpath.bedGraph;
        # -pc = calculates coverage of intervals from left point of a pair reads to the right point
        # exit with error if bedGraph wasn't made
		[ -e $fpath.bedGraph ] || exit 1;

        # 5. Convert to bigWig
        sort -k1,1 -k2,2n $fpath.bedGraph > $fpath.sorted;
        bedGraphToBigWig $fpath.sorted $CHRSZ $fpath.bw;
        
        [ -e $fpath.bw ] || exit 1;
        rm $fpath.sorted;
    fi
done

exit 0;

# if missing, make directory for peak call output
[ -d $FASTQDIR/macs2 ] || mkdir $FASTQDIR/macs2;

## 6. Call peaks
# Specify control(s) to be used
Controls=$FASTQDIR/*WCE*_R1.fastq.gz
for c in $Controls
do
	cpath=${c%_R1.fastq*};
	for i in ${FASTQDIR}/*_R1.fastq.gz
	do
		fpath=${i%_R1.fastq*};
		[ $fpath == $cpath ] && continue;
		outfn=${fpath}_${cpath##*/};

		# -f = input is BAM format
		# -g = genome size. hs=homo sapiens
		# -n = name of outputs
		# -B = output bedGraphs
		# -q = q-value cutoff for peaks
		# --keep-dup = keep all "duplicate" reads in a given location
		macs2 callpeak -t $fpath.bam -c $cpath.bam -f BAMPE -g hs -B -q 0.01 \
		--keep-dup all -n $outfn --outdir $FASTQDIR/macs2;
	done
done
