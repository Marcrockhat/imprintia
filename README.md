# Imprintia

Imprintia is a pipeline designed to discover imprinted genes from DNA-seq and RNA-seq data.
It is flexible in choosing the imprinting threshold and is designed to be compatible for diploid species or species with different ploidy level. It requires parental DNA-seq data and offspring RNA-seq data for imprinting detection.

Dummy files are included for testing the pipeline and as examples of input format.

## Dependencies

This pipeline requires these packages in your Linux system path:

	R 3.0 or above
	Bedtools
	Samtools
	Java 11
	igvtools
	python 3.0 or above (3.5 is installed by default in Ubuntu 14)
	vcf-to-tab
	HTSeq

# Quick Installation

If you're in hurry and don't care what Java 11 will do to your system:

	sudo apt-get install oracle-java11-installer oracle-java11-set-default

Other requisities installation:

	sudo apt-get install r-base bedtools samtools vcftools
	
For HTSeq:

	sudo apt-get install build-essential python2.7-dev python-numpy python-matplotlib python-pysam python-htseq
	
For Freebayes, you will need to compile it by yourself. Please follow their instruction in their page (https://github.com/ekg/freebayes):

	git clone --recursive git://github.com/ekg/freebayes.git

# Quick Guide 

Please tweak this configuration header in the imprintia.sh and adjust it to your system:

	#Configuration: please replace PATH with your sofware path
	#Example
	#samdex=samtools
	#covbed=./genomeCoverageBed
	#freebay=~/freebayes/bin/freebayes
	#vctab=vcf-to-tab
	#igvcount=igvtools

usage:

	./imprintia.sh Abamfile1.bam Bbamfile2.bam AxB.bam BxA.bam reference.fasta annotation.gff3

The first two parameters are bam alignment files to call the parental genotypes from parent A (female) and parent B (male).
The second two parameters are the crosses from both parents.
Reference genome with fasta format is needed. Bias can be introduced if reference used is more similar to one of the parents. A hybrid reference from both parents or from unrelated accession is recommended.
The last parameter are annotation file (GFF3 format) where region of interest will be analyzed for imprinting. This annotation can be changed by any region the user is interested to be analyzed.

You can also type:

	./imprintia.sh -h
	
For future help for the input requirement.

Imprintia.sh will call all the script located in the same folder. For manual running of each script, please see below.

# For Manual Step by Step Use

## SNP Calling from parental genome

Example in bash script:

This step will produce coverage file for each parents.

	samtools view -b A.bam | genomeCoverageBed -ibam crub.chrom.sizes -d > Acov.csv

or this step too, will produce coverage file for each parents.

	genomeCoverageBed -ibam A.bam -g crub.chrom.sizes -d > Acov.txt

This step is needed to call the SNPs from one parental genome:

	freebayes -i -X -u --min-coverage 20 --min-alternate-total 10 -q 30 -b A.bam -v Asnp.vcf -f reference.fasta

then this step is to produce a list of SNPs in easy to read tabulated format.

	cat Asnp.vcf | vcf-to-tab > Atab.vcf


## After SNP calling from parental genome was done, run snpmine.R

Please check and tweak the paramater in this script if you are working with a species with different ploidy level.

	make_option(c("-p", "--ploidy"), type="integer", default = 2, help = "input ploidy level (default = 2)", metavar = "integer")
)

For example, change the "default = 2" into default = 4 if you're working with tetraploid species. 

	Rscript ~/git/imprintia/snpmine.R -a Atab.vcf -b Btab.vcf -c Asort.txt -d Bsort.txt -e /media/diskb/rocky/cruaraproj/SNP/exons.gff -y Acomsnp.csv -z Bcomsnp.csv

## Convert the snp call to query file for temporary reference call by invoking:

The query file contains all SNPs with 'enough' evidence to be informative identifier of each parental genotype based.
No non-informative SNPs will not be included in this list for example: RR --, RR AR, A1A1 A1A2, vice versa and so on.

	#tail here is to grap all the line except the header"
	tail -n +2 Acomsnp.csv | awk '{print $2 ":" $3 "-" $3}' > query.txt

## Call the RNA count. First, after the SNP number extraction using  bashrun.sh:

This process is to quantify RNA reads based on their parental origins.

Command example:

	nohup ./bashrun.sh query2.txt cr48x75.rep1.sorted.rmdup.rg.bam crub.chrom.sizes Cr48xCr75r2.csv &

Notes: the genomic reference fasta file must be at the same place and the chrom.sizes file must have ".fai" end extension.
The script below will require igvtools.jar from IGV and need to be called using java. The script will loop through each coordinate in the query files and count how many RNA reads match the parental origin at the specific position. Since this process is direct and no filter is applied in the quantifying process, a good quality RNA alignment files is needed and filtering need to be done in the alignment process.

	#!/bin/bash
	TEST=$(cat "$1")
	for i in $TEST; do
		java -Xmx1024m  -Djava.awt.headless=true -jar ~/IGVTools/igvtools.jar count -w 1 --bases --query $i $2 count.wig $3
		BASES=$(sed -n '4p' count.wig)
		OUTTEXT=$i$'\t'$BASES
		echo $OUTTEXT >> $4
	done

## Invoke:  fixedsnptest.R

This process is needed to create a useable calculation table needed for the next process.

	Rscript ~/git/imprintia/fixedsnptest.R -i ~/git/imprintia/rna/AxB.bam.csv -s ~/git/imprintia/SNP/Asnp.csv -o ~/git/imprintia/rna/AxBcalc.csv

Afterwards, continue with:  stattest.R

## Then: filteringimp.R

This step is used to test imprinting consistency in two RNA-seq input library. For example, one gene can be found to be maternally imprinted in AxB.bam library but have not enough evidence of similar imprinting in the BxA.bam.

You will need to call HTSeq and Deseq.

Before continuing, a comparison table need to be created first by executing:

	python -m HTSeq.scripts.count -t gene -i ID -f bam <alignment_file> <gff_file>

or 

	for i in ./*.bam; do
		python -m HTSeq.scripts.count -f bam -t gene -i ID $i genes.gff > "${i/%.bam}.csv";
		echo $i is finish!;
	done

## And finally: newsummary.R

This final step is to filter out 
Example:

	Rscript /~/git/imprintia/newsummary.R -a summaryAxB.csv -b summaryBxA.csv -c rawAxB.csv -d rawBxA.csv -w whitelist.csv -l alias.csv -o finalreport.csv

whitelist.csv is a whitelist file that provide additional filter for which genes shall pass the imprinting criteria. If you want to assess genomic imprinting in a tissue that surrounded by maternal origin tissue (eg. dicots endosperm tissue), this list should contain genes that are confirmed to be expressed only in the endosperm to avoid maternal tissue contamination bias.

whitelist.csv is a simple list file which contain gene ID, each gene per line (please check the example table below).

Invoking:

	cat whitelist.csv
Give you:

	GeneA
	GeneB
	GeneC
	

