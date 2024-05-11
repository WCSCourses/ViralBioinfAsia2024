# Reference Alignment To Consensus Practical

## Course Details

* [Viral Genomics and Bioinformatics - Asia]([https://coursesandconferences.wellcomeconnectingscience.org/event/genomics-and-clinical-virology-20230611/](https://coursesandconferences.wellcomeconnectingscience.org/event/viral-genomics-and-bioinformatics-asia-20240512/))
* 12th-17th June 2024
* OUCRU, Vietnam
* [https://github.com/WCSCourses/ViralBioinfAsia2024 ](https://github.com/WCSCourses/ViralBioinfAsia2024 )

## Contact Details

[Dr. Richard Orton](https://www.gla.ac.uk/researchinstitutes/iii/staff/richardorton/)  
[Medical Research Council– University of Glasgow Centre for Virus Research](https://www.gla.ac.uk/research/az/cvr/)  
464 Bearsden Road  
Glasgow  
G61 1QH

E-mail: Richard.Orton@glasgow.ac.uk

## Contents

This practical is associated with a VirtualBox image containing the tools and data sets pre-installed, and an associated lecture on Reference Alignment of High-Throughoput Sequencing (HTS) reads to a reference a sequence, and subsequent variant and consensus calling.

* [1: Setup](#1-setup)
	+ [1.1: Basic read statistics](#11-basic-read-statistics)
* [2: Read Alignment](#2-read-alignment)
	+ [2.1: Indexing the reference sequence](#21-indexing-the-reference-sequence)
	+ [2.2: Aligning the reads to the reference](#22-aligning-the-reads-to-the-reference)
	+ [2.3: Converting SAM to BAM](#23-converting-SAM-to-BAM)
	+ [2.4: Basic alignment statistics](#24-basic-alignment-statistics)
* [3: Alignment on your own](#3-alignment-on-your-own)
* [4: Extra data](#4-extra-data)
* [5: Assembly Visualisation and Statistics Practical](#5-assembly-visualisation-and-statistics-practical)
	+ [5.1: Setup](#51-setup)
	+ [5.2: Summary Statistics with weeSAM](#52-summary-statistics-with-weeSAM)
	+ [5.3: Coverage plot on your own](#53-coverage-plot-on-your-own)
	+ [5.4: Visualisation with Tablet](#54-visualisation-with-tablet)
 

# 1: Setup

In this session, we will first be using a set of Illumina paired end reads which were simulated from a dengue virus (DENV) genome; these simulated reads were created using ART (Huang et al., 2012: [10.1093/bioinformatics/btr708](10.1093/bioinformatics/btr708)). The goal now is to align these reads to a reference genome sequence, with an ultimate goal of creating a consensus sequence for mutation anlysis.

To start off, we will need to move into the correct folder:

```
cd ~/RefAlign/Dengue
```

***Command breakdown:***

* **cd** = change directory 

Next, list the contents of the directory so you can see the files we will be working with:

```
ls
```

You should see the FASTQ paired end read files:

**deng\_sim\_R1.fq**  
**deng\_sim\_R2.fq**

And also two FASTA reference sequence files:

**deng1.fasta**  
**deng3.fasta**

We will be aligning the paired end reads to the two reference sequences in turn. The reference sequences represent the 1 and 3 serotypes of DENV, and we will use the alignment results to determine which serotype the sample contains, and to also highlight the importance of selecting an appropriate reference.
 
## 1.1: Basic read statistics

We will first use a tool called prinseq to count the number of reads in each file. As these are paired end reads, there should be one read from each read pair in each file – and hence the same number of reads in each file. We will also use prinseq to output statistics on the read lengths, but prinseq itself can do much much more.

```
prinseq-lite.pl -stats_info -stats_len -fastq deng_sim_R1.fq -fastq2 deng_sim_R2.fq
```

***Command breakdown:***

1.	**prinseq-lite.pl** is the name of the program
2.	**-stats\_info** tells prinseq to output basic stats on the reads (number of reads and bases)
3.	**-stats\_len** tells prinseq to output basic stats on read lengths (min, max, mean etc)
4.	**-fastq deng\_sim\_R1.fq** the name of the 1st FASTQ file
5.	**-fastq2 deng\_sim\_R2.fq** the name of the 2nd FASTQ file in the pair

### Common Issue
* A common issue here is not entering the prinseq command on one line in the terminal - you should only use the enter key at the end of the command to execute it.
* Another common issue is typos - check the command carefully if you get an error - it is likely you have mispelled a file or argument

***
### Questions
**Question 1** – How many reads and bases are in the read files 1 and 2?

**Question 2** – What is the average (mean) length of the reads? 
***

The statistics are split into those for the first FASTQ file of the read pair (e.g. stats\_info, stats\_len, etc) and those for the second FASTQ file of the read pair (e.g. stats\_info2, stats\_len2, etc), and should look a bit like this:

```
stats_info	bases	5325000
stats_info	reads	35500
stats_info2	bases	5325000
stats_info2	reads	35500
stats_len	max	150
stats_len	mean	150.00
stats_len	median	150
stats_len	min	150
stats_len	mode	150
stats_len	modeval	35500
stats_len	range	1
stats_len	stddev	0.00
stats_len2	max	150
stats_len2	mean	150.00
stats_len2	median	150
stats_len2	min	150
stats_len2	mode	150
stats_len2	modeval	35500
stats_len2	range	1
stats_len2	stddev	0.00

```

Paired read files should always have the same number of lines/reads (the ordering of the reads in each file is also critical), so if your two paired files have a different number of reads, something has gone wrong (e.g. filtering/trimming went wrong and corrupted the output, or maybe files from different samples are being used). 
 
# 2: Read Alignment

There are many tools available to align reads onto a reference sequence: bwa, bowtie2, minimap2, bbMap, to name but a few.

We will be using [BWA](http://bio-bwa.sourceforge.net) to align our paired end reads to a reference sequence and output a [SAM (Sequence Alignment Map)](https://samtools.github.io/hts-specs/SAMv1.pdf) file. The SAM file contains the result of each read’s alignment to the given reference sequence. 

## 2.1: Indexing the reference sequence

First, we need to create a BWA index of the reference sequence. Tools such as BWA need to index the sequence first to create a fast lookup (or index) of short sequence seeds within the reference sequence. This enables the tools to rapidly align millions of reads:

```
bwa index deng1.fasta
```

If you list (ls) the contents of the directory, you should see the BWA index files, they will all have the prefix deng1.fasta, and will have extensions such as **.amb**, **.ann**, **.bwt**, **.pac**, and **.sa**.

```
ls
```

## 2.2: Aligning the reads to the reference

Next, we want to align our reads to the reference sequence using the BWA mem algorithm:

```
bwa mem -t 4 deng1.fasta deng_sim_R1.fq deng_sim_R2.fq > denv1.sam
```

***Command breakdown:***

1. **bwa** = the name of the program we are executing
2. **mem** = the BWA algorithm to use (recommended for illumina reads > 70nt)
3. **-t 4** = use 4 computer threads
4. **deng1.fasta** = the name (and location) of the reference genome to align to
5. **deng\_sim\_R1.fq** = the name of read file 1
6. **deng\_sim\_R2.fq** = the name of read file 2
7. **>** = direct the output into a file
8. **denv1.sam** = the name of the output SAM file to create 

Overall, this command will create an output file called 1b.sam in the current directory, which contains the results (in SAM format) of aligning all our reads to the reference sequence deng1.fasta.

When bwa has finished (and your prompt comes back), check that the SAM file has been created.

```
ls
```

There should now be a file called **denv1.sam** in the directory.

### Common issue
A common mistake is not waiting for your previous command to finish, and entering the next command into the terminal before the prompt has returned. You need to wait until the **manager@ViralGenomics24** command prompt returns before entering the next command - the bwa alignment can sometimes take a few minutes.

## 2.3: Converting SAM to BAM

Typically, a SAM file contains a single line for each read in the data set, and this line stores the alignment result of each read (reference name, alignment location, CIGAR string, the read sequence itself, quality, etc).

SAM files are in a text format (which you can open and view if you like: head 1b.sam), but can take up a lot of disk storage space. It is good practice to convert your SAM files to BAM (Binary Alignment Map) files, which are compressed binary versions of the same data, and can be sorted and indexed easily to make searches faster. We will use [samtools](https://samtools.github.io) to convert our SAM to BAM, and sort and index the BAM file:

```
samtools sort denv1.sam -o denv1.bam
```

```
samtools index denv1.bam
```

***Command breakdown:***

1.	The first command tells samtools to **sort** the SAM file, and to also output (**-o**)the sorted data in BAM format to a file called **1b.bam**
2.	We then use samtools to **index** the BAM file 1b.bam (indexing [which relies on sorted data] enables faster searches downstream).


There should now be two new files in the directory called: 

**1b.bam** (the BAM file)  
**1b.bam.bai** (the BAM index file) 

Now let’s list (ls) the contents of the directory to check we have our new files, and also check out their sizes:

```
ls -lh
```

***Command breakdown:***
* **-l** tells the list (**ls**) command to give the output in a long list format, whilst the **h** tells it to provide file sizes in a human readable format, this is the 5th column, which will have the size of each file in a format such as 8.2M (M for megabytes) or 9.5G (G for gigabytes).

***

### Questions
**Question 3** – How big is the SAM file compared to the BAM file?

***

**NB:** If your SAM file is 0B (i.e. 0 bytes, empty) then something went wrong with the bwa alignment step, so restart from there. If you SAM file is fine (i.e. >0), but your BAM file is 0B (i.e. empty), then something went wrong with your SAM to BAM conversion so re-do that step. 

We don’t need our original SAM file anymore (as we have the BAM file now) so we remove (rm) the SAM file 1b.sam:

```
rm denv1.sam
```

## 2.4: Basic alignment statistics

One common thing to check is how many reads have been aligned (or mapped) to the reference, and how many are not aligned (or unmapped). Samtools can report this for us easily, utilising the aligner SAM flags you learnt about in the previous session.

**Reminder:** the 2nd column in the SAM file contains the flag for the read alignment. If the flag includes the number 4 flag in its makeup then the read is unmapped, if it doesn’t include the number 4 in it's makeup then it is mapped.

### Number of unmapped reads
```
samtools view -c -f4 denv1.bam
```

***Command breakdown***

1.	**samtools view** = to view the file 1b.bam
2.	**–c** = count the read alignments
3.	**–f4** = only include read alignments that do have the unmapped flag 4

### Number of mapped read alignments:
```
samtools view -c -F4 denv1.bam
```

***Command breakdown***

1.	**samtools view** = to view the file 1b.bam
2.	**–c** = count the read alignments
3.	**–F4** = skip read alignments that contain the unmapped Flag 4 

***
### Questions

**Question 4** – how many reads are mapped to the deng1.fasta genome?

**Question 5** – how many reads are unmapped?
***

Technically, the above command gives the number of mapped read **alignments** not reads. A read could be mapped equally well to multiple positions (one will be called the primary alignment, and others secondary alignments [sam flag 256]), or a read could be split into two parts (e.g. spliced) with one part being the primary alignment and the others supplementary [sam flag 2048]

So to get the true number of mapped reads you need to count only the alignments that do not have flags 4 (unmapped), 256 (not primary), and 2048 (supplementary) = 4 + 256 + 2048 = 2308

### Number of mapped reads

```
samtools view -c -F4 -F256 -F2048 denv1.bam
```

or summing up the F flag values together:

```
samtools view -c -F2308 denv1.bam
```

For small RNA viruses, secondary and supplementary alignments tend to be rare, but it is important to know the distinction between mapped **reads** and mapped read **alignments**.

# 3: Alignment on your own

You now need to use bwa to align the reads to the deng3.fasta reference sequence – later in the visualisation and summary statistics section we will be comparing the 1a vs 1b alignment results.

You need to work out the commands yourself based on the previous commands for the deng1.fasta reference. 

Here is a reminder of the commands you used for DENV1 which you will need to adapt. 

**NB:** Essentially, you will want to change the reference name in the bwa command, and all of the SAM/BAM filenames in the bwa and samtools commands from denv1 to denv3.

```
bwa index deng3.fasta
```
```
bwa mem -t 4 deng3.fasta deng_sim_R1.fq deng_sim_R2.fq > denv3.sam
```
```
samtools sort denv3.sam -o denv3.bam
```
```
samtools index denv3.bam
```
```
rm denv3.sam
```
```
samtools view -c -f4 denv3.bam
```
```
samtools view -c -F2308 denv3.bam
```

***
### Questions

**Question 6** – how many reads are mapped to the deng3.fasta genome?

**Question 7** – how many reads are unmapped?

**Question 8** – which reference assembly has the most mapped reads: deng1 or deng3? Therefore, which reference sequence is better (1 or 3)?
***

# 4: Extra Data

If you are looking for something extra to do, there are additional data sets located in the folder:

### ~/RefAlign/Extra/

There are two subfolders in this directory: mystery and mystery2

These are mystery samples, combine all the given references sequences in the **refs** subfolder into one file using the “cat” command, align the reads to that combined reference and then determine what the virus in each sample is.
 
# 5: Assembly Visualisation and Statistics Practical

In this practical, we will be checking our reference assembly from the previous session. We will use tools to generate summary statistics of the depth and breadth of the coverage across the genome, coverage plots, and visualisation of our assembly using tools such as Tablet and weeSAM. Later sessions of the course will cover how to call the consensus sequence and variants.

## 5.1: Setup

In the previous session, you should have bwa aligned the paired reads onto two different DENV genomes (serotypes 1 and 3).

This should have resulted in two BAM files in your Dengue folder, lets check:

```
cd ~/RefAlign/Dengue
```

```
ls
```

You should see (amongst others):

**denv1.bam**  
**denv1.bam.bai**  
**denv3.bam**  
**denv3.bam.bai**  

Along with the two reference sequences:

**deng1.fasta**  
**deng3.fasta**  

We need all these files to proceed, so if you don’t have them – ask for help and we can copy across pre-computed versions.

## 5.2: Summary Statistics with weeSAM

We previously used samtools to count the number of mapped and unmapped reads (using samtools view -c commands), which suggested that DENV3 was a better reference sequence for our sample based on a greater number of mapped reads, but let’s explore this is more detail using a tool called weeSAM: https://github.com/centre-for-virus-research/weeSAM

weeSAM analyses a SAM or BAM file, generates a graphical coverage plot, and reports a range of summary statistics such as:

* **Ref_Name**: The identifier of the reference.
* **Ref_Len**: The length in bases of each reference.
* **Mapped\_Reads**: Number of reads mapped to each reference.
* **Breadth**: The number of sites in the genome covered by reads.
* **%\_Covered**: The percent of sites in the genome which have coverage.
* **Min\_Depth**: Minimum read depth observed.
* **Max\_Depth**: Max read depth observed.
* **Avg\_Depth**: Mean read depth observed.
* **Std\_Dev**: Standard deviation of the mean (Avg_Depth).
* **Above\_0.2_Depth**: Percentage of sites which have greater than 0.2 * Avg_Depth.
* **Above\_1_Depth**: Percentage of sites which are above Avg_Depth.
* **Above\_1.8_Depth**: Percentage of sites which have greater than 1.8 * Avg_Depth.
* **Variation\_Coefficient**: The mean of Std_Dev of the mean.

The Average Depth (Avg_Depth) is perhaps the most important field, along with Breadth which will tell you how much of the genome is covered by aligned reads. But the fields such as Std\_Dev and Above_0.2_Depth can give an indication of the variability in the coverage across the genome.

Let’s run weeSAM on our samples:

```
weeSAM --bam denv1.bam --html denv1
```

An explanation of this command is:

1.	**weeSAM**: the name of the program we are using
2.	**--bam**: flag to signify input bam file
3.	**denv1.bam**: the name of our bam file to analyse
4.	**--html**: flag to signify output html file
5.	**1b**: the name prefix to use for the output

If you list the contents of the directory you should see that a folder called **denv1\_html\_results** has been created:

```
ls
```

Inside this folder is a HTML file that we can view in a web browser (like Firefox or Chrome), the HTML file has the summary statistics and coverage plot so lets take a look and open the html file: 

```
firefox denv1_html_results/denv1.html
```

You should see something like this:

![](https://github.com/WCSCourses/ViralBioinfAsia2024/blob/main/Modules/RefAlign/v2024_denv1_cov.png)


***
### Questions
**Question 9** – what is the average depth of coverage across the deng1 reference genome?
***

Now let’s view the coverage plot by clicking on the hyperlink (blue and underlined) in the Ref_Name column, you should see a coverage plot similar to this:

![](https://github.com/WCSCourses/ViralBioinfAsia2024/blob/main/Modules/RefAlign/v2024_denv3_cov.png)

The x-axis represents the genome position, whilst the y-axis represents the Depth of Coverage at each genome position. 

**NB:** Although our reference file is called deng3.fasta, the actual sequence itself inside the file is called DengueVirus3\_NC\_004175.2 (you can check for yourself if you want to: head –n1 deng3.fasta)

Although you do expect variation in coverage across the genome, the numerous regions of near zero coverage suggest that the DENV1 reference is not ideal, and the aligner has struggled to effectively map reads onto it in this regions – presumably because the reference is too divergent from the viral population in the sample at these regions. 

**Close the weeSAM and Firefox windows before proceeding!**

### Common issue
A common issue here is due to the fact that we have launched firefox from the terminal (wihtout running it background - see advanced linux commands). In order to get our command prompt back (the **manager@ViralGenomics24**) we need to close the firefox window down, the prompt should then return.

## 5.3: Coverage plot on your own

Your task now is to run weeSAM on the deng3.bam file. So you will need to adapt the previous weeSAM command that you used to change the input BAM and output file names

***
### Questions
**Question 10** – what is the average depth of coverage across the deng3 reference genome?

**Question 11** – how does the coverage plot of deng1 compare to deng3? Do you think it is better?
***

## 5.4. Visualisation with Tablet

[Tablet](https://ics.hutton.ac.uk/tablet/) is a tool for the visualisation of next generation sequence assemblies and alignments. It goes beyond simple coverage plots, and allows you to scroll across the genome, zoom into errors of interests, highlight mutations to the reference, and investigate the assembly.

Tablet requires three files:

1.	A bam file, e.g. denv3.bam
2.	A bam index file, e.g. denv3.bam.bai
3.	A reference sequence file: e.g. deng3.fasta

To launch Tablet, type:

```
tablet
```

**NB:** You will not be able to use this command line for other things until you have closed down tablet – but you can open another command line window (or tab) if you want to leave tablet open and do other things.

You should see the Tablet graphical user interface:

![](https://github.com/WCSCourses/GCV23/blob/main/modules/ReferenceAlignment/Tablet.png)

**NB:** Sometimes a small popup window also appears, giving information on how to correctly cite Tablet, with a brief countdown timer.

We want to load in our read alignment from the DENV1 genome. So **Click** on the **Open Assembly** button on the top menu bar.

![](https://github.com/WCSCourses/GCV23/blob/main/modules/ReferenceAlignment/Tablet2.png)

This will launch the Open Assembly window, **Click** **Browse** and then **navigate** to your **~/RefAlign/Dengue** folder and **Select** the **denv.bam** file for **Primary Assembly**. Afterward, **Click** **Browse** and **select** the **deng3.fasta** file for **Reference/Consensus File**, before **Clicking** **Open**.

![](https://github.com/WCSCourses/GCV23/blob/main/modules/ReferenceAlignment/Tablet3.png)

After loading you should see the message **-select a contig to begin visualisation-** along with a list of contigs in the left hand panel. In our analysis, we have used a single sequence (the DENV3 reference sequence), so our contig list only has one entry (the contig DengueVirus3\_NC\_004175.2), **click** on this entry.

![](https://github.com/WCSCourses/GCV23/blob/main/modules/ReferenceAlignment/Tablet4.png)

**NB:** Although our reference file is called deng3.fasta, the actual sequence itself inside the file is called DengueVirus3\_NC\_004175.2 (you can check for yourself if you want to: head –n1 deng3.fasta)

**NB:** We only have one contig as our reference sequence only consisted of one sequence (the DENV3 genome). However, you can align reads to a reference containing multiple sequences, such as the human genome consisting of multiple separate chromosome sequences, or a segmented virus such as influenza consisting of multiple separate segment sequences, or all of the contigs generated from a metagenomics data set. 

Tablet should now load the entire BAM file for visualisation. You can use the **scrollbars** to move across the genome and view all the reads aligned against the reference.

![](https://github.com/WCSCourses/GCV23/blob/main/modules/ReferenceAlignment/Tablet5.png)

### Read Display
In the read display, As, Cs, GS and Ts are represented with different colour blocks, Variants are highlighted with Red Text and a different shading, Deletions are represented with Red Asterisks, whilst the location of Insertions is highlighted with red boxes.

**NB:** Like insertions, soft clipping at the end of the reads are also highlighted with red boxes. Soft clipping is where the aligner has decided to discount a portion of the read (at the read’s beginning or end) to improve the alignment.

You can easily jump about the BAM alignment by **clicking** within the **Coverage Overview** window, and the read display will be updated to show this region.

### Variants
One of the (many) useful features of Tablet is the ability to highlight variants. **Slide** the **Variants Slider** all the way to the right hand side to highlight variants. If you now scroll along the genome, you should be able to easily spot consensus level mutations (as virtually every read at a position will have a mutation) and also spot minority variants.

![](https://github.com/WCSCourses/GCV23/blob/main/modules/ReferenceAlignment/Tablet6.png)

**NB:** Minority variants could be real viral mutations from the viral population or be errors introduced by RT-PCR or the sequencer itself.

***
### Questions

**Question 12:** Can you find a genome position that has a consensus level mutation?
Hint: hold the mouse over a mutation and the genome location will be reported above the read display in red text
***
 
### Colours Schemes

Tablet also has a few other colour schemes for visualisation, accessed through the “Colour Schemes” tab at the top. Try a few out, perhaps the most commonly used schemes are:

1.	Nucleotide: this is the default one: As, Cs, Gs, and Ts represented with different colours
2.	Direction: reads aligned in the forward direction are highlighted in light blue, whilst those in the reverse direction are highlighted in dark blue
3.	Variants: represents As, Cs, Gs, and Ts with grey, and highlights any mutations with red.

### Exit

Remember that you need to close Tablet down in order to get your command line back.

Either click on the red cross in the top left hand corner, or click the Tablet icon (red circle) (located above Open Assembly) and select Exit Tablet.



 




