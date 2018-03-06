# imprintia

Imprintia is a pipeline designed to discover imprinted genes from DNA-seq and RNA-seq data.
It featured felxibility in determining imprinting threshold and are designed to be compatible for diploid or other species with different ploidy level. It requires parental DNA-seq data and offspring RNA-seq data for imprinting detection.

# Quick Guide 

# SNP Calling from parental genome

Example in bash script

samtools view -b Cr48sort.bam | genomeCoverageBed -ibam crub.chrom.sizes -d > 48cov.csv
genomeCoverageBed -ibam Cr48sortstamp.bam -g crub.chrom.sizes -d > 48cov.txt
samtools view -b Cr75sort.bam | genomeCoverageBed -ibam crub.chrom.sizes -d > 75cov.csv
genomeCoverageBed -ibam Cr75sortstamp.bam -g crub.chrom.sizes -d > 75cov.txt
cat Cr48snp.vcf | vcf-to-tab > cr48filtab.vcf
cat Cr75snp.vcf | vcf-to-tab > cr75filtab.vcf

After SNP calling from parental genome was done, run snpmine.R

Rscript ~/git/imprintia/snpmine.R -a Co1719tab.vcf -b Co1979tab.vcf -c Co1719_11sort.txt -d Co1979_9sort.txt -e /media/diskb/rocky/cruaraproj/SNP/crubellaexons.gff -y Co1719comsnp.csv -z Co1979comsnp.csv

Convert the snp call to query file for temporary reference call by invoking:

#tail here is to grap all the line except the header
tail -n +2 compfil48snp.csv | awk '{print $2 ":" $3 "-" $3}' > query.txt

Call the RNA count. First, after the SNP number extraction using  bashrun.sh:
Command example: nohup ./bashrun.sh query2.txt cr48x75.rep1.sorted.rmdup.rg.bam crub.chrom.sizes Cr48xCr75r2.csv &
Notes: the genomic reference fasta file must be at the same place and the chrom.sizes file must have ".fai" end extension.

#!/bin/bash
TEST=$(cat "$1")
for i in $TEST; do
	java -Xmx1024m  -Djava.awt.headless=true -jar ~/IGVTools/igvtools.jar count -w 1 --bases --query $i $2 count.wig $3
	BASES=$(sed -n '4p' count.wig)
	OUTTEXT=$i$'\t'$BASES
	echo $OUTTEXT >> $4
done

Invoke:  fixedsnptest.R
Rscript ~/git/imprintia/fixedsnptest.R -i ~/git/imprintia/rna/Co1719/rep1/AxB.bam.csv -s ~/git/imprintia/SNP/Asnp.csv -o ~/git/imprintia/rna/Co1719/rep1/AxBcalc.csv

Afterwards, continue with:  stattest.R

Then: filteringimp.R

HTSeq and Deseq

Before continuing, this table need to be created first by executing:

python -m HTSeq.scripts.count -t gene -i ID -f bam <alignment_file> <gff_file>

or 

for i in ./*.bam; do
	python -m HTSeq.scripts.count -f bam -t gene -i ID $i crubgenes.gff > "${i/%.bam}.csv";
	echo $i is finish!;
done

Afterwards, continue with contaminationcheck.R

And finally: newsummary.R
Rscript /~/git/imprintia/newsummary.R -a summaryAxB.csv -b summaryBxA.csv -c rawAxB.csv -d rawBxA.csv -w whitelist.csv -l alias.csv -o finalreport.csv

whitelist.csv should contain gene ID, each gene per line.
