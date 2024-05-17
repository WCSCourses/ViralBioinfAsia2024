# Step 1

#cp ../DeNovo/SARS-CoV-2/SRR21065613_1.fastq \
#../DeNovo/SARS-CoV-2/SRR21065613_2.fastq \
#../DeNovo/SARS-CoV-2/Ref/SARS-CoV-2.fa \
#.
basequal=15
baselen=60
read1=$(basename $1)
read2=$(basename $2)
read3=$(basename $3)

# Cleaning the Reads
trim_galore -q $basequal --length $baselen --paired $read1 $read2
 
# Running Prinseq
prinseq-lite.pl -fastq $read1 -fastq2 $read2 -out_format 3 -min_len 60 -min_qual_score 15

# Index the Reference
echo "indexing"
bwa index $read3 

# Map the Reads and Creating Sam and BAM file 
echo "Mapping"
bwa mem -t 4 $read3 $read1 $read2> SARS_Group1.sam

samtools sort SARS_Group1.sam -o SARS_Group1.bam
samtools index SARS_Group1.bam
ls -lhr
rm SARS_Group1.sam

# Create Mapped Reads 
samtools view -c -F4 SARS_Group1.bam

# Generate the Consensus Sequence
samtools mpileup -aa -A -d 0 -Q 0 SARS_Group1.bam | ivar consensus -p SARS_Group1_consensus -t 0.4


#Look for variants
lofreq call -f $read3 -o SARS_CoV.vcf SARS_Group1.bam

conda activate snpeff
snpEff -ud 0 $read3 SARS_CoV.vcf > SARS_CoV_snpeff.vcf
conda deactivate

