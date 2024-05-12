
# Metagenomics

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

This practical is associated with a VirtualBox image containing the tools and data sets pre-installed, and an associated lecture on Metagenomics encompassing two parts: kmer based methods on reads and blast based methods on contigs.

* [1: Kmer based metagenomics with Kraken2](#1-kmer-based-metagenomics-with-kraken2) 
* [2: Run kraken2](#2-run-kraken2)
* [3: Kraken2 on your own](#3-kraken2-on-your-own)
* [4:Contig metagenomics with diamond](4-contig-metagenomics-with-diamond)
* [5:Diamond on your own](5-diamond-on-your-own)

# 1: Kmer based metagenomics with Kraken2

In this session, we will be working with two Illumina metagenomic data sets that have already been published and are available for download on the [NCBI Short Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra), there are two samples:

* Human: location = ~/Metagenomics/Human/
* Vampire bat: ocation = ~/Metagenomics/Vampire/

We will be using a tool called [Kraken2](https://ccb.jhu.edu/software/kraken2/) to analyse the paired end FASTQ reads for each sample. Kraken2 is the newest version of Kraken, a taxonomic classification system using exact k-mer matches to achieve high accuracy and fast classification speeds. This classifier matches each k-mer within a query sequence (i.e. a read or contig) to the lowest common ancestor (LCA) of all genomes within the database containing the given k-mer.

We will be broadly following the published Kraken metagenomic protocol:

**Metagenome analysis using the Kraken software suite**  
Lu et al. (2022). Nature Protocols volume 17, pages 2815–2839  
[https://www.nature.com/articles/s41596-022-00738-y](https://www.nature.com/articles/s41596-022-00738-y)  

Normally we would be using the "standard" kraken database - which can be downloaded from the [Kraken database download page](https://benlangmead.github.io/aws-indexes/k2). The standard database includes RefSeq archaea, bacteria, viruses, plasmid complete genomes, UniVec Core and the most recent human reference genome, GRCh38. 

However, due to space (and computational power) constraints on the Virtual Machine we will be using the Viral RefSeq database, which is very small and fast, but as the dtaabase is built only using viral sequences it does come at a cost of false positives. In addition, kraken is nucleotide based, and the database typically only contains with genome sequence for each viral species, so may miss very divergent viruses.

To set things up, first change directory (cd) to directory where we will be working:

```
cd ~/Metagenomics/Human/
```

Next, list the contents of the directory so you can see the files we will be working with:

```
ls
```

You should see the FASTQ paired-end read files which are gzipped (compressed):

**SRR533978\_1.fq.gz**  
**SRR533978\_2.fq.gz**

# 2: Run kraken2


```
kraken2 --db ~/kraken2_database/ --threads 6 --minimum-hit-groups 3 --report-minimizer-data --paired SRR533978_1.fastq.gz SRR533978_2.fastq.gz --output kraken_output.txt --report kraken_report.txt
```

***Command breakdown:***

1. **kraken2** = the name of the program we are executing
2. **--db ~/kraken2_database/** = the location of the kraken2 database we are using 
3. **--threads 6** = use 6 computer threads
4.  **--report-minimizer-data** = forces Kraken 2 to provide unique k-mer counts per classification
5. **--minimum-hit-groups 3** = for increased classification precision (i.e., fewer false positives).
6. **--paired SRR533978\_1.fastq.gz SRR533978\_2.fastq.gz** = the name of the paired input FASTQ files
7. **--output kraken_output.txt** = the name of the kraken output file to create
8. **--report kraken_report.txt** = the name of the kraken report output file to create

**NB:** The --minimum-hit-groups flag specifies the minimum number of ‘hit groups’ needed to make a classification call. Hit groups are overlapping k-mers sharing the same minimizer. Kraken 2 uses minimizers to compress the input genomic sequences, thereby reducing storage memory needed and run time. In this example, we increase the minimum number of hit groups from the default two groups to three groups for increased accuracy. See [Kraken2 manual](https://github.com/DerrickWood/kraken2/wiki/Manual)

Overall, this command will output two files - our kraekn\_output.txt file and our kraken\_report.txt file. The kraken\_output.txt is large and not really human readable whereas the kraken_report.txt is a human readable summary that is tab-delimited, with one line per taxon. The fields of the output, from left-to-right, are as follows:

1. Percentage of reads covered by the clade rooted at this taxon
2. Number of reads covered by the clade rooted at this taxon
3. Number of reads assigned directly to this taxon
4. Number of minimizers in read data associated with this taxon (new)
5. An estimate of the number of distinct minimizers in read data associated with this taxon (new)
6. A rank code, indicating (U)nclassified, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies. All other ranks are simply '-'.
7. NCBI taxonomy ID
8. indented scientific name

You can explore the report file manually using a command like **more**, esepcially if you understand the taxonomy. However, a more visual way of viewing the results is via a Krona plot:

```
ktImportTaxonomy -q 2 -t 3 -s 4 kraken_output.txt -o kraken_krona.html -tax ~/Diamond/
```
***Command breakdown:***

1. **ktImportTaxonomy** = the name of the program we are executing
2. **-q 2** = column 2 of input files to use as query ID.
3. **-t 3** = column 3 of input files to use as taxonomy ID
4. **-s 4** = column 4 of input files to use as score. 
5. **kraken_output.txt** = the name of the input file - whcih is the kraken output file
6. **-o kraken_krona.html** = the name of the Krona html file to create
7. **-tax /db/kronatools/taxonomy** = the location of the NCBI taxonomy files for Krona

This will create an new html output file which you should see in your directory:

```
ls
```

We can open this file with Firefox:

```
firefox kraken_krona.html
```

***
### Questions
**Question 1** – What viruses have the highest read counts in the sample?
***

**NB:** Alternatively, the html file can be downloaded via the MobaXterm file browser on the left hand side of the window onto your local machine and opened there

## 2.1 What is a krona plot?

Krona is an interactive visualization tool for exploring the composition of metagenomes within a Web browser. A Krona plot is a html file that can be opened by a web browser (such as Firefox, Chrome, Safari and Internet Explorer/Edge). Krona uses multilevel pie charts to visualize both the most abundant organisms and their most specific classifications (Fig. 1). Rather than hiding lower ranks in its overview, Krona hides low-abundance organisms, which can be expanded interactively. 

**Krona: Interactive Metagenomic Visualization in a Web Browser**  
Ondov et al. (2013)  
DOI 10.1007/978-1-4614-6418-1_802-1  
[https://link.springer.com/content/pdf/10.1007/978-1-4614-6418-1_802-1.pdf](https://link.springer.com/content/pdf/10.1007/978-1-4614-6418-1_802-1.pdf)

An example Krona plot is shown below:

![](https://github.com/rjorton/OIE2023/blob/main/krona.png)

Some of the key features of a Krona plot:

1. You can double click on any taxon (i.e. any slice of the pie) to view only that taxon and it's children
2. Clicking/Selecting a taxon will display how many sequences (i.e. reads or contigs) have been assigned to that taxon in the top right corner under "count" (this is actually a hyperlink to a file containing all the sequence IDs assigned to this taxon). The taxonomy ID is also displayed as a hyperlink to the corresponding NCBI page.
3. The defauly colouring of the krona plot is set to distinguish taxons apart - however the colouring can be changed based on the 'score' (for BLAST results this corresponds to e-value) by selecting the "Color by Avg. log e-value" - this would highlight the taxons with the highest score


## 2.2: Bas congo virus / Tibrovirus congo

The Human sample came from the following paper where deep sequencing was used to discover a novel rhabdovirus (Bas-Congo virus, or BASV) associated with a 2009 outbreak of 3 human cases of acute hemorrhagic fever in Mangala village, Democratic Republic of Congo (DRC), Africa.

**A Novel Rhabdovirus Associated with Acute Hemorrhagic Fever in Central Africa**  
Grard et al. 2012  
PLoS Pathog. 2012 Sep; 8(9): e1002924.  
[https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3460624/](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3460624/)

The viral species is now called 'Tibrovirus congo': [https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=1987017](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=1987017)

# 3: Kraken2 on your own

If you have time, try processing the second sample the vampire bat and reusing and adapting the commands for the human sample - although there is no need to align the human genome first.

The sample is one sample from the following paper:

**Using noninvasive metagenomics to characterize viral communities from wildlife**  
Molecular Ecology Resources Volume19, Issue 1, January 2019, Pages 128-143  
Bergner et al. 2018  
[https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.12946](https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.12946)

# 4: Contig metagenomics with diamond

First lets move into the correct folder:

```
cd ~/Metagenomics/Human/
```

If we have time we can run spades to de novo assemble the reads into large conrigs:

```
metaspades.py -1 SRR533978_1.fastq.gz -2 SRR533978_2.fastq.gz -t 6 -o ./spades
```

Otherwise we can use the contigs I already assembled when I ran spades before the course, located here:

```
~/Metagenomics/Human/Pre/contigs.fasta
```

Now we run diamond blastx on the contigs, on a viral refseq database (for speed and limited resources):

```
diamond blastx --threads 6 --outfmt 6 --db ~/Diamond/viral-refseq-prot --out contigs_diamond.txt -q contigs.fasta 
```
Now we generate a krona plot of the results:

```
ktImportBLAST contigs_diamond.txt -o contigs_diamond.html -tax ~/Diamond/

firefox contigs_diamond.html
```

# 5:Diamond on your own

```
cd ~/Metagenomics/VIZIONS
```

Pre assembled contigs here:

```
~/Metagenomics/VIZIONS/Pre/contigs.fasta
```

