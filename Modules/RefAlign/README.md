# Reference Alignment To Consensus Practical

## Course Details

* [Viral Genomics and Bioinformatics - Asia](https://coursesandconferences.wellcomeconnectingscience.org/event/viral-genomics-and-bioinformatics-asia-20240512/)
* 12th-17th May 2024
* OUCRU, Vietnam
* [https://github.com/WCSCourses/ViralBioinfAsia2024 ](https://github.com/WCSCourses/ViralBioinfAsia2024)

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
* [6: Consensus sequence generation](#6-consensus-sequence-generation)
* [7: Variant calling](#7-variant-calling)
* [8: Group practical](#8-group-practical)

# 1: Setup

In this session, we will first be using a set of Illumina paired end reads which were simulated from a dengue virus (DENV) genome; these simulated reads were created using ART (Huang et al., 2012: [10.1093/bioinformatics/btr708](10.1093/bioinformatics/btr708)). The goal now is to align these reads to a reference genome sequence, with an ultimate goal of creating a consensus sequence for mutation anlysis.

To start off, we will need to move into the correct folder:

```
cd ~/RefAlign/Dengue
```

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

We will be aligning the paired end reads to the two reference sequences in turn. The reference sequences represent serotypes 1 and 3 of DENV, and we will use the alignment results to determine which serotype the sample contains, and to also highlight the importance of selecting an appropriate reference.
 
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

As we are using simulated reads they are very uniform (all the same length) and also very good quality. For this practical we are skipping read trimming and QC as the reads are such good quality and as read QC is covered in a separate practical. But, in normal circumstances you should always do QC on your reads!
 
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

Overall, this command will create an output file called denv1.sam in the current directory, which contains the results (in SAM format) of aligning all our reads to the reference sequence deng1.fasta.

When bwa has finished (and your prompt comes back), check that the SAM file has been created.

```
ls
```

There should now be a file called **denv1.sam** in the directory.

### Common issue
A common mistake is not waiting for your previous command to finish, and entering the next command into the terminal before the prompt has returned. You need to wait until the **manager@ViralGenomics24** command prompt returns before entering the next command - the bwa alignment can sometimes take a few minutes.

## 2.3: Converting SAM to BAM

Typically, a SAM file contains a single line for each read in the data set, and this line stores the alignment result of each read (reference name, alignment location, CIGAR string, the read sequence itself, quality, etc).

SAM files are in a text format (which you can open and view if you like: head denv1.sam), but can take up a lot of disk storage space. It is good practice to convert your SAM files to BAM (Binary Alignment Map) files, which are compressed binary versions of the same data, and can be sorted and indexed easily to make searches faster. We will use [samtools](https://samtools.github.io) to convert our SAM to BAM, and sort and index the BAM file:

```
samtools sort denv1.sam -o denv1.bam
```

```
samtools index denv1.bam
```

***Command breakdown:***

1.	The first command tells samtools to **sort** the SAM file, and to also output (**-o**)the sorted data in BAM format to a file called **denv1.bam**
2.	We then use samtools to **index** the BAM file denv1.bam (indexing [which relies on sorted data] enables faster searches downstream).


There should now be two new files in the directory called: 

**denv1.bam** (the BAM file)  
**denv1.bam.bai** (the BAM index file) 

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

We don’t need our original SAM file anymore (as we have the BAM file now) so we remove (rm) the SAM file denv1.sam:

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

1.	**samtools view** = to view the file denv1.bam
2.	**–c** = count the read alignments
3.	**–f4** = only include read alignments that do have the unmapped flag 4

### Number of mapped read alignments:
```
samtools view -c -F4 denv1.bam
```

***Command breakdown***

1.	**samtools view** = to view the file denv1.bam
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

You now need to use bwa to align the reads to the deng3.fasta reference sequence – later in the visualisation and summary statistics section we will be comparing the denv1 vs denv3 alignment results.

You need to work out the commands yourself based on the previous commands for the deng1.fasta reference. 

**NB:** Essentially, you will want to change the reference name in the bwa command, and all of the SAM/BAM filenames in the bwa and samtools commands from denv1 to denv3.

Here is a reminder of the commands you used for DENV1 **which you will need to adapt**. 

```
bwa index deng1.fasta
```
```
bwa mem -t 4 deng1.fasta deng_sim_R1.fq deng_sim_R2.fq > denv1.sam
```
```
samtools sort denv1.sam -o denv1.bam
```
```
samtools index denv1.bam
```
```
rm denv1.sam
```
```
samtools view -c -f4 denv1.bam
```
```
samtools view -c -F2308 denv1.bam
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

These are mystery samples, combine all the given references sequences in the **refs** subfolder into one file using the “cat” command, align the reads to that combined reference and then determine what the virus in each sample is - based on the number of reads mapping to each reference sequence. An easy command to use here is:

```
samtools idxstats INPUT.bam
```
**NB:** you need to change the name INPUT to whatever your BAM file is actually called

This will report the number of mapped read alignments (alignments not reads), the reference sequence length, and number of unmapped reads on each of the sequences in the reference file used for the alignment e.g:

```
Ref1	Seqlen 12345	0
Ref2	Seqlen 654	12
Ref3	Seqlen 2	0
*	Seqlen 65276	0
```

idxstats represents unmapped reads in two forms. First reads where BOTH members of the pair are unapped are recorded on the last line where the Reference name is '*'. Secondly, reads that are unmapped but whose pair did map to a reference sequence are recorded in the 3rd column of the corresponding reference sequences. For example, in the example above, Ref2 had 654 reads that did map, and 12 reads that were techincally unmapped but those 12 reads each had their pair (from the paired end reads) map to the designated reference sequence.

Another useful function of samtools (and there are lots) is the flagstat command:

```
samtools flagstat INPUT.bam
```

This will report the number of read alignments, number of primary, secondary and supplementary alignments and much more.
 
***
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
5.	**denv1**: the name prefix to use for the output

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

![](https://github.com/WCSCourses/ViralBioinfAsia2024/blob/main/Modules/RefAlign/v2024_denv1_cov2.png)

The x-axis represents the genome position, whilst the y-axis represents the Depth of Coverage at each genome position. 

**NB:** Although our reference file is called deng1.fasta, the actual sequence itself inside the file is called DengueVirus1\_NC\_004177.1 (you can check for yourself if you want to: head –n1 deng1.fasta)

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

![](https://github.com/WCSCourses/ViralBioinfAsia2024/blob/main/Modules/RefAlign/Tablet.png)

**NB:** Sometimes a small popup window also appears, giving information on how to correctly cite Tablet, with a brief countdown timer.

We want to load in our read alignment from the DENV3 genome. So **Click** on the **Open Assembly** button on the top menu bar.

![](https://github.com/WCSCourses/ViralBioinfAsia2024/blob/main/Modules/RefAlign/Tablet2.png)

This will launch the Open Assembly window, **Click** **Browse** and then **navigate** to your **~/RefAlign/Dengue** folder and **Select** the **denv3.bam** file for **Primary Assembly**. Afterward, **Click** **Browse** and **select** the **deng3.fasta** file for **Reference/Consensus File**, before **Clicking** **Open**.

![](https://github.com/WCSCourses/ViralBioinfAsia2024/blob/main/Modules/RefAlign/v2024_tablet_select.png)

After loading you should see the message **-select a contig to begin visualisation-** along with a list of contigs in the left hand panel. In our analysis, we have used a single sequence (the DENV3 reference sequence), so our contig list only has one entry (the contig NC\_004175.2), **click** on this entry.

![](https://github.com/WCSCourses/ViralBioinfAsia2024/blob/main/Modules/RefAlign/Tablet4.png)

**NB:** We only have one contig as our reference sequence only consisted of one sequence (the DENV3 genome). However, you can align reads to a reference containing multiple sequences, such as the human genome consisting of multiple separate chromosome sequences, or a segmented virus such as influenza consisting of multiple separate segment sequences, or all of the contigs generated from a metagenomics data set. 

Tablet should now load the entire BAM file for visualisation. You can use the **scrollbars** to move across the genome and view all the reads aligned against the reference.

![](https://github.com/WCSCourses/ViralBioinfAsia2024/blob/main/Modules/RefAlign/Table5.png)

### Read Display
In the read display, As, Cs, GS and Ts are represented with different colour blocks, Variants are highlighted with Red Text and a different shading, Deletions are represented with Red Asterisks, whilst the location of Insertions is highlighted with red boxes.

**NB:** Like insertions, soft clipping at the end of the reads are also highlighted with red boxes. Soft clipping is where the aligner has decided to discount a portion of the read (at the read’s beginning or end) to improve the alignment.

You can easily jump about the BAM alignment by **clicking** within the **Coverage Overview** window, and the read display will be updated to show this region.

### Variants
One of the (many) useful features of Tablet is the ability to highlight variants. **Slide** the **Variants Slider** all the way to the right hand side to highlight variants. If you now scroll along the genome, you should be able to easily spot consensus level mutations (as virtually every read at a position will have a mutation) and also spot minority variants.

![](https://github.com/WCSCourses/ViralBioinfAsia2024/blob/main/Modules/RefAlign/v2024_tablet_mutations.png)

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


# 6: Consensus sequence generation

We have now aigned our DENV sample to the dengue serotype 3 reference genome sequence, and now we want to call a consensus sequence.

What is a consensus sequence? At each genome position in the SAM/BAM alignment file, call the most frequent nucleotide (or insertion/deletion) observed in all of the reads aligned at the position. 

In this practical, we will use a tool called [iVar](https://andersen-lab.github.io/ivar/html/manualpage.html) to call the consensus sequence, which utilises the [mpileup](http://www.htslib.org/doc/samtools-mpileup.html) function of samtools.

First, lets make sure we are in the correct directory

```
cd ~/RefAlign/Dengue
```

And now call the consenus for the sample using iVar:

```
samtools mpileup -aa -A -d 0 -Q 0 denv3.bam | ivar consensus -p denv3_consensus -t 0.4
```

Breaking this command down, there are two parts:

1. samtools [mpileup](http://www.htslib.org/doc/samtools-mpileup.html) which essentially outputs the base and indel counts for each genome position
	* **-aa** = output data for absolutely all positions (even zero coverage ones)
	* **-A** = count orphan reads (reads whose pair did not map)
	* **-d 0** = override the maximum depth (default is 8000 which is typically too low for viruses)
	* **-Q 0** = minimum base quality, 0 essentially means all the data
2. ivar [consensus](https://andersen-lab.github.io/ivar/html/manualpage.html) - this calls the consensus - the output of the samtools mpileup command is piped '|' directly into ivar
	* -p denv3_consensus = prefix with which to name the output file
	* -t 0.4 = the minimum frequency threshold that a base must match to be used in calling the consensus base at a position. In this case, an ambiguity code will be used if more than one base is > 40% (0.4). See [iVar manual](https://andersen-lab.github.io/ivar/html/manualpage.html)

By default, iVar consensus uses a minimum depth (-m) of 10 and a minimum base quality (-q) of 20 to call the consensus; these defaults can be changed by using the appropriate arguments. If a genome position has a depth less than the minimum, an 'N' base will be used in the consensus sequence by default.

iVar will output some basic statistics to the screen such as:

```
#DO NOT ENTER THIS - IT IS AN EXAMPLE OF AN IVAR OUTPUT:
[mpileup] 1 samples in 1 input files
[mpileup] Max depth set to maximum value (2147483647)
Minimum Quality: 20
Threshold: 0.4
Minimum depth: 10
Minimum Insert Threshold: 0.8
Regions with depth less than minimum depth covered by: N
Reference length: 10707
Positions with 0 depth: 0
Positions with depth below 10: 5
```

and when it has finished (and your prompt returns) you should see our consensus sequence (denv3_consensus.fa) in the directory:

```
ls
```

which you can view the sequence via the command line (we will be covering variants later):

```
more denv3_consensus.fa 
```

***
### Questions

**Question 13** - are there any Ns in the generated consensus sequence at the start? why do you think that is?
***

**Question 14** - try copying and pasting the created consensus sequence into [NCBI Blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi) - what is the closest sample on GenBank?
***


# 7: Variant calling

Viruses, and in particular RNA viruses, can exist as complex populations consisting of numerous variants present at a spectrum of frequencies – the so called viral quasispecies. Although we have created a consensus sequence (which typically considers mutations at a frequency >50% in the sample) using iVar, we do not yet know anything about the mutations within the sample. Furthermore, it is often necessary to go beyond the consensus, and investigate the spectrum of low frequency mutations present in the sample.

iVar itself could be used to call variants (using the [iVar variants](https://andersen-lab.github.io/ivar/html/manualpage.html) command). But here we will be using a slightly more advanced variant caller called [LoFreq](https://github.com/CSB5/lofreq) to call the low (and high) frequency variants present in the sample BAM file. LoFreq uses numerous statistical methods and tests to attempt to distinguish true low frequency viral variants from sequence errors. It requires a sample BAM file and corresponding reference sequence that it was aligned to, and creates a [VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf) file as an output.

First, lets make sure we are in the correct folder to work on Sample1:

```
cd ~/RefAlign/Dengue
```
To use LoFreq enter this command:

```
lofreq call -f deng3.fasta -o denv3.vcf denv3.bam
```

Breaking this command down:

* **lofreq**: the name of the program we are using
* **call**: the name of the function within LoFreq we are using – call variants
* **-f deng3.fasta**: the reference file name and location (path)
* **-o deng3.vcf**: the output VCF file name to create
* **deng3**.bam: the input BAM file name

Now lets open the VCF file created by LoFreq:

```
more denv3.vcf
```

The outputted VCF file consists of the following fields:

* CHROM: the chromosome – in this case the DENV3 sequence NC_001475.2 
* POS: the position on the chromosome the variant is at
* ID: a ‘.’ but LoFreq can be run with a database of known variants to annotate this field
* REF: the reference base at this position
* ALT: the alternate base (the mutated base) at this position
* QUAL: LoFreq’s quality score for the variant
* FILTER: whether it passed LoFreq’s filters e.g. PASS
* INFO: Detailed Information
	* DP=1248; depth = 1248
	* AF=0.995192; Alt Frequency (Mutation Frequency) = 99.5192%
	* SB=0; Strand Bias test p-value
	* DP4=0,1,604,638: Coverage of the ref base in Fwd and Rev, and the alt base in Fwd and Rev

***
### Questions

**Question 15** – how many consenus level (i.e AF > 0.5) are there in the sample? Are there any subconsenus mutations (i.e. AF < 0.5) 
***

LoFreq simply calls the reference position and mutation present, it does not characterise the effect of the mutation in terms of whether the mutation is synonymous or nonsynonymous etc. To do that we will use a program called [SnpEff](https://github.com/pcingola/SnpEff) which is run on LoFreq’s outputted VCF file and creates a new annotated VCF file.

First we need to activate a conda environment snpeff is installed in. Due to a conflict on the VM, snpEff was installed in its own special environment, to use it we need to activate it:

```
conda activate snpeff
```
Then run snpeff on our vcf file to annotate it:

```
snpEff -ud 0 NC_001475.2 denv3.vcf > denv3_snpeff.vcf
```

Breaking this command down:

* **snpEff**: the name of the program we are using
* **-ud 0**: Set upstream downstream interval length to 0
* **NC_001475.2**: the reference file name
* **deng3.vcf**: the input vcf file name
* **deng3_snpeff.vcf**: the annotated output vcf file name

Setting -ud 0 stops SnpEff from characterising mutations located near (but not within) a gene as being located in their UTRs. Typically viral genomes are compact and genes are separated by few bases, in such cases a mutation in one gene could also be characterised as being in the UTR region of a neighbouring gene (as SnpEff was initially built for human analyses) – try running SnpEff without the -ud 0 and compare the results if you want, you should see multiple annotations for each mutation.

Now we can view the annotated vcf file created by SnpEff:

```
more denv3_snpeff.vcf
```

The mutations will now have annotations added at the end of the Info field (the last field, the 8th field), e.g in Sample1 you should see this nonsynonymous (missense) mutation at genome position 644 which corresponds to a Leu to Val amino acid mutation at codon 184 (Leu184Val) in the polyprotein (POLY) on DENV3:

```
DO NOT ENTER THIS - IT IS AN EXAMPLE OF A MUTATION IN THE VCF!!!
NC_001475.2	644	.	C	G	34996.0	PASS DP=1012;AF=0.994071;SB=6;DP4=0,3,496,510;ANN=G|missense_variant|MODERATE|POLY|DV3_gp1|transcript|DV3_gp1|protein_coding|1/1|c.550C>G|p.Leu184Va
l|550/10173|550/10173|184/3390||

```

The C (4th column) to G (5th column) mutation at genome position 644 (2nd column) corresponds to position 550 (out of 10173) within the polyprotein (POLY) gene which corresponds to codon 184 (out of 3390) within POLY.

**NB:** SnpEff includes many pre-built databases – for many viruses you may need to build the SnpEff database first by downloading and processing a GenBank file, see the documentation [here](https://pcingola.github.io/SnpEff/)

When we are finished using snpEff we need to deactivate the conda environment:

```
conda deactivate
```

# 8: Group practical

For the group practical, we want you to work in groups on a single computer and analyse one (or both) of the below dengue virus samples. Tasks:

* Align the reads to the designated reference sequence using bwa
* Convert the SAM to BAM etc
* Count the number of mapped and unmapped reads
* Create a coverage plot
* What is the average depth and breadth of the sample
* Create a consensus sequences with ivar - BLAST it on GenBank - what is the closest sequence?
* Use LoFreq to call mutations and create a vcf file
* Use snpEff to characterise the mutations - find one synonymous and one non-synonymous mutation


**Dengue-S**
* Virus: Dengue virus - serotype 2
* Country: Seychelles
* Location: ~/RefAlign/Dengue-S
* Reference file: dengue_refseq.fasta
* FASTQ files: SRR10230796_1.fastq  and SRR10230796_2.fastq
* Paper Title: [Complete Genome Sequences of Dengue Virus Type 2 Epidemic Strains from Reunion Island and the Seychelles](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6979306/)

**Dengue-B**
* Virus: Dengue virus - serotype 2
* Country: Bangladesh
* Location: ~/RefAlign/Dengue-B
* Reference file: dengue_refseq.fasta
* FASTQ files: SRR12901070_1.fastq  and SRR12901070_2.fastq
* Paper Title: [Genome Sequence of a Dengue Virus Serotype 2 Strain Identified during the 2019 Outbreak in Bangladesh](https://journals.asm.org/doi/10.1128/mra.01246-20)

One Friday, we can try and write a bash script to process thee automatically.

# 9: Moving forward

Typically when aligning to a good reference sequence, you are essentially just re-using the same commands again and again. Index the reference, align the reads to the reference, convert the SAM to BAM, index the BAM, count the number of reads, etc etc. This is where the power of BASH scripting and bioinformatics workflows comes into force. All you are doing is changing the names of the input and output files, the commands themselves are not changing.

```
bwa index REF.fasta
```
```
bwa mem -t 4 REF.fasta Read1.fq Read2.fq > NAME.sam
```
```
samtools sort NAME.sam -o NAME.bam
```
```
samtools index NAME.bam
```
```
rm NAME.sam
```
```
samtools view -c -f4 NAME.bam
```
```
samtools view -c -F2308 NAME.bam
```

One very good bioinformatics workflow for processing viral data is [ViralRecon](https://github.com/nf-core/viralrecon) which is written in nextflow. It is capable of automatically processing your FASTQ reads to cerate consensus sequences, you need to install nextflow and ViralRecon and then prepare input files specifying the input file names and locations.

Another alternative is the Galaxy web servers for processing HTS data: e.g. [https://usegalaxy.eu](https://usegalaxy.eu)




