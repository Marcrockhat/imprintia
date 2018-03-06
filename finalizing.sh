#!/bin/bash
##  MRH 31032016
##  finalizing.sh v0, script to finalize output of Imprintia if two replicates were used
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
$1 : replicate 1 compiled result
$2 : replicate 2 compiled
$3 : first cross bam
$4 : second cross bam
$5 : reference genome
$6 : annotation file
"

Rscript newsummary.R -
