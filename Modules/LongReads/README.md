# Long read data analysis

## Course Details

* [Viral Genomics and Bioinformatics - Asia](https://coursesandconferences.wellcomeconnectingscience.org/event/viral-genomics-and-bioinformatics-asia-20240512/)
* 12th-17th May 2024
* OUCRU, Vietnam
* [https://github.com/WCSCourses/ViralBioinfAsia2024 ](https://github.com/WCSCourses/ViralBioinfAsia2024)

## Contact

[Richard Orton](https://www.gla.ac.uk/schools/infectionimmunity/staff/richardorton/)   
[MRC-University of Glasgow Centre for Virus Research](https://www.gla.ac.uk/research/az/cvr/)  
464 Bearsden Road  
Glasgow  
G61 1QH  
UK  
E-mail: Richard.Orton@glasgow.ac.uk  

## Contents

* [1: HCMV Nanopore](#1-hcmv-nanopore)
	+ [1.1: HCMV Nanopore on your own](#11-hcmv-nanopore-on-your-own)
* [2: ARTIC Nanopore SARS-CoV-2 Reference Alignment](#2-artic-nanopore-sars-cov-2-reference-alignment)
 	+ [2.1: Setup and Data](#21-setup-and-data)
	+ [2.2: ARTIC consensus sequence generation walkthrough](22-ARTIC-consensus-sequence-generation-walkthrough)
	+ [2.3: Generating ARTIC consensus sequences yourself](23-generating-artic-consensus-sequences-yourself)


Previously in the course you would of learnt about the FASTQ format, what SAM and BAM files, and been aligning primarily illumina reads to reference sequences to call a consensus sequence. This session will build upon this, tweaking the steps to adapt them for handling nanopore reads, and also for handling the large number of overlapping amplicons used in the ARTIC protocols.

## 1: HCMV Nanopore

Our first nanopore data set is from two different Human Cytomegalovirus (HCMV) samples, one from urine and one from the lung but from different patients. We will do a walk-through analysis of the urine sample, and then you will need to apply the commands to the lung sample.

One of the best features of nanopore reads is the lengths the reads can reach. So first, let's run the same prinseq command from the illumina session to get an idea of the average and maximum read lengths

```
prinseq-lite.pl -stats_info -stats_len -fastq SRR23882357_urine.fastq 
```

**Question:** how many reads are there?

**Question:** what is the average read length?

**Question:** what is the longest read length?

We can use the tool [NanoPlot](https://github.com/wdecoster/NanoPlot) to further explore read lengths (through various plots) but also explore the quality of the reads:

```
NanoPlot -t 6 --fastq SRR23882357_urine.fastq --only-report
```

Then we can open the NanoPlot report file we just created:

```
firefox NanoPlot-report.html 
```

**Question:** what is the average quality of reads? how does this compare to yesterday's illumina data?

The next step is to align the reads (**SRR23882357\_urine.fastq**) to the reference sequence (**urine.fasta**). Here, we need to use an aligner that can handle the higher error rate of nanopore reads - such as [minimap2](https://github.com/lh3/minimap2):

```
minimap2 -x map-ont -t 4 urine.fasta SRR23882357_urine.fastq  -a -o urine.sam
```
**-x map-ont** tells minimap2 it is nanopore data  
**-t 4** means use 4 computer threads  
**-a** means output the data in SAM format  
**-o** is the output file name  

Now we have created a SAM file, we can do our now standard SAM to BAM conversion, sort, index, and creation of a coverage plot:

```
samtools sort urine.sam -o urine.bam

samtools index urine.bam

weeSAM --bam urine.bam --html urine

firefox urine_html_results/urine.html 
```

**Question:** how many read are mapped? what is the average coverage across the genome?

As we now have a BAM file, we could still use iVar to generate a consensus sequence like we did for the illumina data yesterday, i.e:

```
samtools mpileup -aa -A -d 0 -Q 0 urine.bam | ivar consensus -q 10 -p urine_consensus -t 0.4
```

However, as nanopore is suspectible to sequencing errors (particularly indels) at homopolymers it is better to use a consensus caller capable of correcting these sequencing errors such as [medaka](https://github.com/nanoporetech/medaka):

```
medaka_consensus -i SRR23882357_urine.fastq -d urine.fasta -m r941_min_high_g360 -t 6 
```
**-i** is our original fastq files
**-d** is our reference sequence 
**-m** is the error correction model to use which encompasses the nanopore chemistry (r941), hardware (minion rather than promethion), base calling mode (high) & guppy version (g360)

This will create a folder call **medaka** with a consensus sequence file called **consensus.fasta** which we can view with:

```
more medaka/consensus.fasta
```

You could try BLASTing the sequence to find the closest match (this may take a while on NCBI)


If you want (**YOU DO NOT NEED TO DO THIS**), you could try aligning the consensus sequence to the original sequence and exaning difference in AliView:

```
cat urine.fasta medaka/consensus.fasta > all_seqs.fa

mafft all_seqs.fa > aligned.fa

aliview aligned.fa
```

## 1.1: HCMV Nanopore on your own

Now try working the same steps for the lung sample:

**Reference: lung.fasta**  
**Reads: SRR23882358_lung.fastq**

## 2: ARTIC Nanopore SARS-CoV-2 Reference Alignment

This nanopore tutorial is based on the [ARTIC](https://artic.network) network's [nCoV-2019 bioinformatics protocol](https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html) which we will use to create consensus genome sequences for a number of MinION SARS-CoV-2 samples. 

## 2.1: Setup and Data

The 'artic-ncov2019' [Conda](https://docs.conda.io/en/latest/) environment has already been installed on the [VirtualBox](https://www.virtualbox.org) Ubuntu virtual machine (VM) that you are using for this course, but we need to 'activate' the artic-ncov2019 Conda environment each time we want to use the ARTIC pipeline:

```
conda activate artic-ncov2019
```

Next, lets change directory (```cd```) into the folder where the example nanopore data is located for this practical:

```
cd ~/LongReads/SARS-COV-2/

```

If you list the contents of this directory you should see a number of files and folders:

```
ls
```

You should see the following barcode folders (06, 07 and 12) representing the 3 samples on the run that we will be analysing:

* **barcode06**
* **barcode07**
* **barcode12**


Typically the FASTQ data for each sample is stored in multiple files of around 4000 reads each. For barcode06, you should see 25 different FASTQ files, numerically labelled at the end of their filename from 0 to 24:

```
ls barcode06
```


## 2.2: ARTIC consensus sequence generation walkthrough

The first sample we will be working with is barcode06. The ARTIC bioinformatics protocol has two distinct steps:

1. **artic guppyplex** - combines all a samples FASTQ reads into a single file and size filters them (it can also perform a quality score check, which is not needed here as the reads are already split into pass and fail folders based on quality by guppy)
2. **artic minion** - aligns the reads to the [Wuhan-Hu-1](https://www.ncbi.nlm.nih.gov/nuccore/MN908947) reference sequence, trims the amplicon primer sequences from the aligned reads, downsamples amplicons to reduce the data, and creates a consensus sequence utilising [medaka](https://github.com/nanoporetech/medaka) for variant calling to correct for common MinION errors (such as those associated with homopolymer regions).

First we will create a folder to work in and store our output files:

```
mkdir Results
```

Then we will move into the folder to work:

```
cd Results
```

Now we will run artic guppyplex on sample barcode06:

```
artic guppyplex --skip-quality-check --min-length 400 --max-length 700 --directory ../barcode06 --prefix cvr124a
```

Breaking this command down:

* **artic gupplyplex** = the name of the program/module to use (installed as part of conda environment)
* **--skip-quality-check** = don't filter reads based on quality (our reads are already filtered)
* **--min-length 400** = minimum read length to accept is 400 bases
* **--max-length 700** = maximum read length to accept is 700 bases
* **--directory** = PATH to input directory containing FASTQ reads to process
* **--prefix** = output name prefix to label output file (I choose cvr124a to signify batch 124a from the CVR)

This should create an output FASTQ file called **cvr124a_barcode06.fastq**:

```
ls
```

**artic minion** - next we will run artic minion using the FASTQ file created above:


```
artic minion --normalise 200 --threads 4 --scheme-directory ~/artic-ncov2019/primer_schemes --read-file cvr124a_barcode06.fastq --medaka --medaka-model r941_min_high_g360 nCoV-2019/V2 barcode06
```

Breaking this command down:

* **artic minion** = the name of the program/module to use (installed as part of conda environment)
* **--normalise 200** = normalise (downsample) each amplicon so there are only 200 reads in each direction (forward and reverse) - this enables the variant and consensus calling to complete relatively quickly
* **--threads 4** = the number of computer threads to use (depends on how powerful your machine is, the more the better)
* **--scheme-directory** = path to the artic primer scheme directory (installed on the VM as part of the conda environment)
* **--read-file** = the name of the input FASTQ file to align
* **--medaka** = use the medaka variant caller (instead of nanopolish)
* **--medaka-model** = the error correction model - nanopore chemistry, hardware, base calling mode & guppy version
* **nCoV-2019/V2** = the primer scheme to use for amplicon primer trimming - this folder is located in the scheme\_directory so on this VM this corresponds to the folder ~/artic-ncov2019/primer_schemes/nCoV-2019/V2
* **barcode06** = the output prefix name to label output files, this can be anything you want such as the sample name or barcode number or anything you want

Overall, this artic minion command uses the aligner [minimap2](https://github.com/lh3/minimap2) to align the reads to the [Wuhan-Hu-1](https://www.ncbi.nlm.nih.gov/nuccore/MN908947) reference sequence, [samtools](http://www.htslib.org) to convert and sort the SAM file into BAM format, custom [artic](https://github.com/artic-network/fieldbioinformatics) scripts for amplicon primer trimming and normalisation (downsampling), [medaka](https://github.com/nanoporetech/medaka) for variant calling, and custom [artic](https://github.com/artic-network/fieldbioinformatics) scripts for creating the consensus sequence using the reference and VCF files, and then masking low coverage regions with Ns. This will create the following key files (amongst many others), all starting with the prefix **barcode06** in this instance:

* **barcode06.sorted.bam** - BAM file containing all the reads aligned to the reference sequence (there is no amplicon primer trimming in this file)
* **barcode06.trimmed.rg.sorted.bam** - BAM file containing normalised (downsampled) reads with amplicon primers left on - this is the file used for variant calling
* **barcode06.primertrimmed.rg.sorted.bam** - BAM file containing normalised (downsampled) reads with amplicon primers trimmed off
* **barcode06.pass.vcf.gz** - detected variants that PASSed the filters in VCF format (gzipped)
* **barcode06.fail.vcf** - detected variants that FAILed the filters in VCF format
* **barcode06.consensus.fasta** - the consensus sequence of the sample

We can view the FASTA consensus sequence via the command line:

```
more barcode06.consensus.fasta
```

You could trt BLASTing it or uploading it to [USHER](https://genome.ucsc.edu/cgi-bin/hgPhyloPlace) to see what the consensus is similar to.

Alternatively you could aslo try creating a coverage plot with weeSAM - but you will need to deactivate the artic conda environment first

## 2.3: Generating ARTIC consensus sequences yourself

Your task now is to adapt the above artic guppyplex and minion commands to run on samples barcode07 and/or barcode12.


At the end of this session we should deactivate our conda environment:

```
conda deactivate
```









