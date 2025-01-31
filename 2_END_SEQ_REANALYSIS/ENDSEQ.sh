#!/bin/bash

CONDA_PATH=$(conda info | grep -i 'base environment' | awk '{print $4}')
source $CONDA_PATH/etc/profile.d/conda.sh

basenames=$(find ./ -maxdepth 1 -type f -name "*.fastq.gz" | while read F; do basename "$F" | rev | cut -c 13- | rev; done | sort | uniq)
echo $basenames
sleep 10

conda activate BASIC_ALIGNMENT
for i in $basenames; do echo "Quality filtering file: $i"
	fastp --in1 "$i"_R1.fastq.gz --in2 "$i"_R2.fastq.gz --out1 "$i"_tR1.fastq.gz --out2 "$i"_tR2.fastq.gz -Q --dedup --dup_calc_accuracy 4 --thread 8 --html ./FASTQ/"$i".fastp.html --json ./FASTQ/"$i".fastp.json 
	mv "$i"_R1.fastq.gz "$i"_R2.fastq.gz ./FASTQ/; done
	
for i in $basenames; do echo "Aligning $i to reference genome"
	bowtie2 --very-sensitive-local --threads=8 -x /mnt/c/Users/Herman/Desktop/NGS/_BOWTIE_INDEXES/ECOLI/MG1655_GENOME_MASKED_TER -1 "$i"_tR1.fastq.gz -2 "$i"_tR2.fastq.gz -S "$i".sam; done

for i in $basenames; do echo "Creating sorted BAMs for $i"
	samtools view -b -@ 8 --min-MQ 30 -F 0x4 "$i".sam | samtools sort -o "$i".bam -m 1000000000 -@ 8 -
	samtools index "$i".bam -@ 8; done

conda activate DEEPTOOLS
for i in $basenames; do echo "Creating bigwigs for $i"
	bamCoverage \
	-b "$i".bam \
	-o "$i"_for.bigwig \
	-of bigwig \
	--normalizeUsing CPM \
	-bs 1000 \
	--filterRNAstrand reverse \
	--smoothLength 1 \
	--Offset 1 \
	-p max/2 \
	--verbose
	bamCoverage \
	-b "$i".bam \
	-o "$i"_rev.bigwig \
	-of bigwig \
	--normalizeUsing CPM \
	-bs 1000 \
	--filterRNAstrand forward \
	--smoothLength 1 \
	--Offset 1 \
	-p max/2 \
	--verbose
done

for i in $(find ./ -maxdepth 1 -type f -name "*.bigwig" | while read F; do basename $F | rev | cut -c 14- | rev; done | sort | uniq)
	do multiBigwigSummary \
	bins \
	--bwfiles "$i"_1_for.bigwig "$i"_2_for.bigwig \
	-out "$i".npz \
	--binSize 1000 \
	-p 8 \
	--verbose \
	--outRawCounts "$i"_for.bedgraph
	rm "$i"_1_for.bigwig "$i"_2_for.bigwig "$i".npz
	
	multiBigwigSummary \
	bins \
	--bwfiles "$i"_1_rev.bigwig "$i"_2_rev.bigwig \
	-out "$i".npz \
	--binSize 1000 \
	-p 8 \
	--verbose \
	--outRawCounts "$i"_rev.bedgraph
	rm "$i"_1_rev.bigwig "$i"_2_rev.bigwig "$i".npz; done

for i in $(find ./ -maxdepth 1 -type f -name "*.bedgraph" | while read F; do basename $F | rev | cut -c 10- | rev; done | sort | uniq)
	do temp_file="temp.tsv"
	tail -n +2 "$i.bedgraph" > "trimmed_input.tsv"
	awk -F'\t' 'BEGIN { OFS="\t" } { sum=0; for (i=4; i<=NF; i++) sum+=$i; print sum/(NF-3) }' "trimmed_input.tsv" > "$temp_file"
	cut -f 1-3 "trimmed_input.tsv" | paste -d '\t' - "$temp_file" > "$i"_temp.bedgraph
	sort -k2,2 -n "$i"_temp.bedgraph > "$i"_avg.bedgraph
	rm "trimmed_input.tsv" "$temp_file" "$i"_temp.bedgraph; done

conda activate BASIC_ALIGNMENT
for i in $(find ./ -maxdepth 1 -type f -name "*_avg.bedgraph" | while read F; do basename $F | rev | cut -c 14- | rev; done | sort | uniq)
	do bedGraphToBigWig \
	"$i"_avg.bedgraph \
	chrom.size \
	"$i".bigwig; done






