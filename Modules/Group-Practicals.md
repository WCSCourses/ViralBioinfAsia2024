## Ref mapping

- Ref mapping script
	### Essential: 
	- Input: ref.fa file_1.fq (file_2.fq) and output name

	- Processing: 
		- Clean the rears
		- Index the reference
		- Map the reads
	- Output: file.sam(bam)
 
	  
	### Desired: 
	- Check the input
	- Create a sam(bam) with mapped reads only
	- Mapping statistics
	- Consensus sequence
	- Look for variants


	### Advanced:
	- Give the help message
	- Check whether the input files are in the right format
	- Give the option to choose the mapping program (bwa, bowtie2, tanoti, minimap2)
	- Create the sam(bam) file with the same name as the input filename (with an appropriate file extension)
	- Check whether the reads are Illumina/ONT and run the analysis accordingly

