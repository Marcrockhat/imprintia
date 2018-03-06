#!/bin/bash
##  MRH 31032016
##  imprintia.sh v0, script to run all imprintia pipeline
##  This source code file is part of the analysis pipeline Imprintia
##	Copyright (C) 2016 Marcelinus Rocky Hatorangan
## 
##  This program is free software: you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation, either version 3 of the License, or
##  (at your option) any later version.
## 
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
## 
##  You should have received a copy of the GNU General Public License
##  along with this program.  If not, see <http://www.gnu.org/licenses/>.
##
##  Please cite this pipeline if you use it in your work
##
##  Marcelinus Rocky Hatorangan <marcelinusrocky@aol.com>
##  31 March 2016
##
##  Tested on Ubuntu 14.04
##
##  input files:
##  bam files for both parental line (bam 1 and bam 2), reference genome (fasta), chromosome file size, 
##
##  output files:
##  compiled SNP list for both parents

echo "For help type:

imprintia.sh -h

Imprintia
pipeline to identify imprinted genes
usage:
	./imprintia.sh Abamfile1.bam Bbamfile2.bam AxB.bam BxA.bam reference.fasta annotation.gff3
where:
	-h  show this help text"

#Configuration: please fill in the braket in the $() with your sofware path
#samdex=$(samtools)
#covbed=$(genomeCoverageBed)
#freebay=$(~/FreeBayes/bin/freebayes)
#vctab=$(vcf-to-tab)

#Example usage samtools view -b CgA.bam | genomeCoverageBed -ibam crub.chrom.sizes -d > CgAcov.csv
# This "${1/%\.bam/cov}" means replace ".bam" with "cov"
echo "
$1 : parental bam 1
$2 : parental bam 2
$3 : first cross bam
$4 : second cross bam
$5 : reference genome
$6 : annotation file
"

#mkdir ./temp

filename1=$(basename $1)
filename2=$(basename $2)
filename3=$(basename $3)
filename4=$(basename $4)
filename5=$(basename $5)

#Preparing input
echo "$5"
#samtools faidx "$5"
#mv "$5".fai ./temp/"$filename5".size

#Calling SNPs
#genomeCoverageBed -ibam "$1" -g ./temp/"$filename5".size -d > ./temp/"${filename1/%\.bam/cov}".csv
#genomeCoverageBed -ibam "$2" -g ./temp/"$filename5".size -d > ./temp/"${filename2/%\.bam/cov}".csv

#~/FreeBayes/bin/freebayes -i -X -u --min-coverage 20 --min-alternate-total 10 -q 30 -b "$1" -v ./temp/"${filename1/%\.bam/snp}".vcf -f "$5"
#~/FreeBayes/bin/freebayes -i -X -u --min-coverage 20 --min-alternate-total 10 -q 30 -b "$2" -v ./temp/"${filename2/%\.bam/snp}".vcf -f "$5"

#cat ./temp/"${filename1/%\.bam/snp}".vcf | vcf-to-tab > ./temp/"${filename1/%\.bam/filtab}".csv
#cat ./temp/"${filename2/%\.bam/snp}".vcf | vcf-to-tab > ./temp/"${filename2/%\.bam/filtab}".csv

#Compiling SNPs
#Rscript snpmine.R -a ./temp/"${filename1/%\.bam/filtab}".csv -b ./temp/"${filename2/%\.bam/filtab}".csv -c ./temp/"${filename1/%\.bam/cov}".csv -d ./temp/"${filename2/%\.bam/cov}".csv -e "$6" -y ./temp/"${filename1/%\.bam/comsnp}".csv -z ./temp/"${filename2/%\.bam/comsnp}.csv" 

#Preparing the file for nucleotide calling
tail -n +2 ./temp/"${filename1/%\.bam/comsnp}".csv | awk '{print $2 ":" $3 "-" $3}' > ./temp/query.txt
./nucleocount.sh ./temp/query.txt "$3" ./temp/"$filename5".size "${filename3}".csv
tail -n +2 ./temp/"${filename2/%\.bam/comsnp}".csv | awk '{print $2 ":" $3 "-" $3}' > ./temp/query.txt
./nucleocount.sh ./temp/query.txt "$4" ./temp/"$filename5".size "${filename4}".csv

#Compiling counts

Rscript fixedsnptest.R -i ./temp/"${filename3}".csv -s ./temp/"${filename1/%\.bam/comsnp}".csv -o ./temp/"${filename3}"calc.csv
Rscript fixedsnptest.R -i ./temp/"${filename4}".csv -s ./temp/"${filename2/%\.bam/comsnp}".csv -o ./temp/"${filename4}"calc.csv

Rscript stattest.R -i ./temp/"${filename3}"calc.csv -o ./temp/"${filename3}"det.csv
Rscript stattest.R -i ./temp/"${filename4}"calc.csv -o ./temp/"${filename4}"det.csv

#Final result
mkdir ./result

Rscript filteringimp.R -a ./temp/"${filename3}"det.csv -b ./temp/"${filename4}"det.csv -c ./result/compiled.csv -r ./result/raw.csv -o ./result/summary.csv
