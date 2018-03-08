# Imprintia

Imprintia is a pipeline designed to discover imprinted genes from DNA-seq and RNA-seq data.
It is flexible in choosing the imprinting threshold and is designed to be compatible for diploid species or species with different ploidy level. It requires parental DNA-seq data and offspring RNA-seq data for imprinting detection.

## Dependencies

This pipeline requires these packages in your Linux system path:

	R 3.0 or above
	Bedtools
	Samtools
	igvtools
	python 3.0 or above
	vcf-to-tab
	HTSeq

# Quick Guide 

usage:

	./imprintia.sh Abamfile1.bam Bbamfile2.bam AxB.bam BxA.bam reference.fasta annotation.gff3

or

	./imprintia.sh --help
	
For future help.
	
# For Manual Step by Step Use

## SNP Calling from parental genome

Example in bash script

	samtools view -b A.bam | genomeCoverageBed -ibam crub.chrom.sizes -d > Acov.csv

or

	genomeCoverageBed -ibam A.bam -g crub.chrom.sizes -d > Acov.txt

then

	cat Asnp.vcf | vcf-to-tab > Atab.vcf

## After SNP calling from parental genome was done, run snpmine.R

	Rscript ~/git/imprintia/snpmine.R -a Atab.vcf -b Btab.vcf -c Asort.txt -d Bsort.txt -e /media/diskb/rocky/cruaraproj/SNP/exons.gff -y Acomsnp.csv -z Bcomsnp.csv

## Convert the snp call to query file for temporary reference call by invoking:

	#tail here is to grap all the line except the header"
	tail -n +2 Acomsnp.csv | awk '{print $2 ":" $3 "-" $3}' > query.txt

## Call the RNA count. First, after the SNP number extraction using  bashrun.sh:

Command example:

	nohup ./bashrun.sh query2.txt cr48x75.rep1.sorted.rmdup.rg.bam crub.chrom.sizes Cr48xCr75r2.csv &

Notes: the genomic reference fasta file must be at the same place and the chrom.sizes file must have ".fai" end extension.

	#!/bin/bash
	TEST=$(cat "$1")
	for i in $TEST; do
		java -Xmx1024m  -Djava.awt.headless=true -jar ~/IGVTools/igvtools.jar count -w 1 --bases --query $i $2 count.wig $3
		BASES=$(sed -n '4p' count.wig)
		OUTTEXT=$i$'\t'$BASES
		echo $OUTTEXT >> $4
	done

## Invoke:  fixedsnptest.R

	Rscript ~/git/imprintia/fixedsnptest.R -i ~/git/imprintia/rna/AxB.bam.csv -s ~/git/imprintia/SNP/Asnp.csv -o ~/git/imprintia/rna/AxBcalc.csv

Afterwards, continue with:  stattest.R

## Then: filteringimp.R

HTSeq and Deseq

Before continuing, this table need to be created first by executing:

	python -m HTSeq.scripts.count -t gene -i ID -f bam <alignment_file> <gff_file>

or 

	for i in ./*.bam; do
		python -m HTSeq.scripts.count -f bam -t gene -i ID $i genes.gff > "${i/%.bam}.csv";
		echo $i is finish!;
	done

## Afterwards, continue with contaminationcheck.R

## And finally: newsummary.R

Example:

	Rscript /~/git/imprintia/newsummary.R -a summaryAxB.csv -b summaryBxA.csv -c rawAxB.csv -d rawBxA.csv -w whitelist.csv -l alias.csv -o finalreport.csv

whitelist.csv should contain gene ID, each gene per line.
