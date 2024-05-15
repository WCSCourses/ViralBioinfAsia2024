## De novo assembly practicals


### [Sreenu Vattipally](https://www.gla.ac.uk/schools/infectionimmunity/staff/sreenuvattipally/)
* MRC-University of Glasgow Centre for Virus Research
* University of Glasgow, G61 1QH
* E-mail: Sreenu.Vattipally@glasgow.ac.uk

### Course information 

* [Viral Genomics and Bioinformatics - Asia](https://coursesandconferences.wellcomeconnectingscience.org/event/viral-genomics-and-bioinformatics-asia-20240512/)
* 12th-17th May 2024
* OUCRU, Vietnam




In the past decade, next-generation sequencing (NGS) has become an important tool for understanding viral evolution. It also has revolutionised many areas of biological research. NGS platforms can now generate large amounts of sequence data in a relatively short time, making them well-suited for studies of genomes, transcriptomes, and epigenetics.

One of the important applications of NGS is whole-genome sequencing (WGS), which can be used to obtain the complete genome sequence of an organism. However, the primary challenge with WGS is that it is difficult to assemble the genome without a reference. 

De novo assembly is the process of assembling a genome without using any reference. This can be a challenge since it requires correctly assembling millions of short DNA sequences (reads) that cover the entire genome. This approach has been used to discover new viruses, characterise known viruses, and study the evolution of viruses.

There are several steps involved in viral de novo assembly, including, quality control, read trimming, assembly, scaffolding of contigs, and gap filling.



In this practical, we will be using three different assemblers i.e. Spades, IDBA UD and ABYSS to assemble the SARS-CoV-2 genome from Illumina reads. 


Create a directory for the analysis
```
mkdir -p ~/DeNovo/SARS-CoV-2
cd ~/DeNovo/SARS-CoV-2/
```

Data for this practical (SRA ID: SRR21065613) are in your home directory's “Sreenu/FQs” folder. Move them to the current directory and uncompress them.

```
mv ~/Sreenu/FQs/SRR21065613_* .
gunzip SRR21065613_*
```


Do a quality check and trimming using trim galore program.

```
trim_galore -q 20 -length 50 --illumina --paired SRR21065613_1.fastq SRR21065613_2.fastq
```

After running trim galore you should see two new files (`SRR21065613_1_val_1.fq` and `SRR21065613_2_val_2.fq`) in your working directory. 

Now create individual directories for Spades, IDBA UD and Abyss along with Reference and Validation folders. 

```
mkdir Spades Abyss IDBA_UD Ref Validate
```

## De novo assembly using Spades  

First, go to the Spades directory you have created above

```
cd ~/DeNovo/SARS-CoV-2/Spades
```

Run spades assembly program with different k-mer sizes.

```
spades.py -k 21,45,73,101 -1 ../SRR21065613_1_val_1.fq -2 ../SRR21065613_2_val_2.fq -o .
```

Here

- -k : different k-mer sizes    
- -1: First reads file
- -2: Second reads file
- -o: output location. Here we are saving in the current directory(.), we can also specify any name.

The program will take some time depending on the CPU power and available memory to complete the job. If everything goes without errors, the final assembled contigs will be stored in contigs.fasta file. In this file, contigs are saved in a sorted manner, with the longest contig first. Remember, these contigs can be positive or negative strands. We have to check the orientation comparing them with the reference sequence.
    
**TODO: What is the length of the longest contig generated?  

## De novo assembly using IDBA_UD 

Now let us try running IDBA UD. 

```
cd ~/DeNovo/SARS-CoV-2/IDBA_UD/
```

Though it is an assembly program for high throughput reads, IDBA_UD works only with fasta files(not fastq). It will not consider the quality of the reads while assembling. So... we have to convert reads from fastq to fasta and submit them to IDBA_UD.

```
fq2fa --merge ../SRR21065613_1_val_1.fq ../SRR21065613_2_val_2.fq reads.fa
idba_ud --mink=21 --maxk=81 step=10 -r reads.fa -o . --num_threads=4
```

Here

- fq2fa: program to convert reads from fastq to fasta
    
- mink: minimum k-mer size
- maxk: maximum k-mer size
- step: k-mer size increment
- r: input reads
- o: output directory
- num_threads: number of processing threads to run the program

The final results will be stored in contig.fa file in the current directory.

**TODO: Check the length of the longest contig.

## De novo assembly using ABYSS 

Go the Abyss directory

```
cd ~/DeNovo/SARS-CoV-2/Abyss/
```

Run Abyss assembler

```
abyss-pe k=51 n=4 in='../SRR21065613_1_val_1.fq ../SRR21065613_2_val_2.fq' name=Abyss-51 B=2G
```
    
Here
- k: k-mer size
- n: number of threads to use
- in: input reads
- name: output name

As you notice from above, the abyss assembler works with one k-mer at a time. Since we do not know the optimal k-mer size for our data it is advisable to use various k-mer sizes. Use the below `for loop` to run the program with different k-mers.

```
for k in 41 51 61 71 81; 
	do 
	abyss-pe k=$k n=4 in='../SRR21065613_1_val_1.fq ../SRR21065613_2_val_2.fq' name=Abyss-$k B=2G; 
done
```

This will run abyss with k-mer sizes from 41 to 81 with the increment of 10 and save the output in corresponding files. To get the assembly stats, run

```
abyss-fac *-contigs.fa
```

Now combine all the contigs generated by different kmer sizes

```
cat Abyss*-unitigs.fa > Abyss-contigs.fa
```



## Assembly validation

```
cd ~/DeNovo/SARS-CoV-2/Validate/

mv ~/Sreenu/Ref/SARS-CoV-2.fa ../Ref/.
mv ~/Sreenu/GFF/SARS-CoV-2.gff ../Ref/.

```

Let us check whether the contigs belong to SARS-CoV-2, if so where they map and what is the percentage of genome covered by them. For this we will be using the ”nucmer” program to compare the contigs to a reference genome.

```
nucmer ../Ref/SARS-CoV-2.fa ../Spades/contigs.fasta -p Spades
mummerplot -p Spades -t png Spades.delta
eog Spades.png
```

This will open a graph of nucmer mapping of all the spades generated contigs, their mapping location and orientation. Here X-axis is the reference genome location and Y-axis is the contigs mapped positions. If a contig is shown in ascending manner it is mapping to reference sequence positive strand. If it is shown in descending manner it is mapping to the negative strand of the reference sequence.  

**TODO: run the nucmer mapping with Abyss and IDBA UD contigs and check their orientation and coverage. 

Now let us use quast program to validate the assemblies. (below command is a single line command)

```
quast.py  -l "Spades, IDBA_UD, Abyss" ../Spades/contigs.fasta ../IDBA_UD/contig.fa ../Abyss/Abyss-contigs.fa -R ../Ref/SARS-CoV-2.fa --genes ../Ref/SARS-CoV-2.gff --reads1 ../SRR21065613_1_val_1.fq --reads2 ../SRR21065613_2_val_2.fq
```

This will create a directory called `quast_results` and save the results. Open them using firefox since they are in html format.

```
firefox quast_results/latest/report.html
```

**TODO: Can you identify which assembler produced better results? and why? Do you see any misassembled blocks? Which program misassembled the contigs?


## Scaffolding of contigs

In this section, we will be doing reference guided scaffolding of the de novo contigs using `contigsMerger` progreom. Each contig is mapped to the reference to find its location and orientation. If there are any overlapping contigs, they will be merged to create longer(super) contigs.

```
cd ~
git clone https://github.com/vbsreenu/ContigsMerger.git
```

Create a separate directory to save the results

```
mkdir  ~/DeNovo/SARS-CoV-2/GapFilling
cd ~/DeNovo/SARS-CoV-2/GapFilling
```

Combine all the contigs into a file and run contigsMerger program.

```
cat ../Abyss/*-contigs.fa ../Spades/contigs.fasta ../IDBA_UD/idba_ud/contig.fa > all_contigs.fa
```

```
~/ContigsMerger/contigsMerger ../Ref/SARS-CoV-2.fa all_contigs.fa N
```

This will create a super-contig combining all the contigs. Unassembled regions will be filled with ‘N’s.

## Gaps filling

If there are no reads in the dataset that cover a region of the genome, there is no way we can fill those gaps computationally. However, if you think there are reads covering the entire genome, running different assemblers with varying k-mer sizes should generate contigs that will fill all the gaps. 

Another alternative way to fill the gaps is by taking a gap corresponding sequence from a close reference sequence. However, after generating a super-contig, reads should be mapped to it and a consensus sequence should be generated. If there are enough reads in the dataset, the consensus sequence will reflect it.

```
~/ContigsMerger/contigsMerger ../Ref/SARS-CoV-2.fa all_contigs.fa
```

This will fill the gaps with the reference sequence. Use the “combinedContig.fa” sequence as a reference and map read using bwa mem. Check the assembly in “tablet” for any indels or mismatches around gap regions.

---

## Group practical
1. Download the fastq data from [NCBI](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&page_size=10&acc=SRR13638088&display=download)
2. Create a directory called "deNovoPractical" in your home directory.
3. Create individual folders for different de novo assembly results. 
4. Trim the reads (check how you can use single-end data with trim_galore)
5. Run spades, idba_ud and abyss assembly programs. (again, check the parameters for running single-end reads)
6. Check the longest contig from each assembly
7. Run BLAST search to check which genome the contig belongs to.
8. Make a reverse complement of the genome if it is aligning to the minus strand
9. Check if there are any structural variations in the assembly by comparing it with the closest reference genome (you have to get the whole genome fasta file and gene features file, i.e gff)
