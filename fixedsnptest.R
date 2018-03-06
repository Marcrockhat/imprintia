#!/usr/bin/env Rscript
##  MRH 18112013
##  
##  fixedsnptest.R v0, This script is used to identify maternal and paternal SNP and calculate the p-value on proportion test based on several methods.
##  
##  This source code file is part of the analysis pipeline Imprintia
##  Copyright (C) 2016 Marcelinus Rocky Hatorangan
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

rm(list = ls())

args <- commandArgs(TRUE)

#Opt section
library('optparse')
option_list = list(
	make_option(c("-i", "--input"), type = "character", default = NULL, help = "parent input", metavar = "character"),
	make_option(c("-s", "--snptable"), type = "character", default = NULL, help = "parent snps", metavar = "character"),
	make_option(c("-o", "--output"), type = "character", default = NULL, help = "output file", metavar = "character"),
#this is a test feature
    make_option(c("-p", "--ploidy"), type="integer", default = 2, help = "input ploidy level (default = 2)", metavar = "integer")
)

opt_parser = OptionParser(option_list= option_list)
opt = parse_args(opt_parser)
#End of opt section

snpdat <- read.table(opt$input, header= FALSE, sep= " ", fill= TRUE)
names(snpdat) = c("NAMES","POS", "NUMA", "NUMC", "NUMG", "NUMT", "NUMN")

#Revert the column names for technical reason, I'm using ATGC as my conventional way to name column
#flipping this order will change the outcome of the analysis since I'm using fixed index which related to the column name for the bases
#A is 1, T is 3, G is 9, and C is 27. If you want to change the order, please rewrite the index which is used to distinguish the maternal and paternal alleles
snpdat <- snpdat[c("NAMES","POS", "NUMA", "NUMT", "NUMG", "NUMC", "NUMN")]
snpdat[is.na(snpdat)] <- 0

#Extract additional information for each locus including parental allele type and gene name.
snpref <- read.table(opt$snptable, header= TRUE, sep= "\t", fill= TRUE)
snpdat$MAT <- snpref$MAT
snpdat$PAT <- snpref$PAT
genecat <- snpref$genecat
exoncat <- snpref$exoncat

#rm(snpref)
matsnp <- as.matrix(snpdat$MAT)
newmatsnp <- as.matrix(snpdat$MAT)
patsnp <- as.matrix(snpdat$PAT)
newpatsnp <- as.matrix(snpdat$PAT)
nallele <- 4
ploidy <- 2

#greed is a variable container which decides whether shared allele between maternal and paternal will be discarded
#if greed equal to 1, the shared allele will be keeped and will be added as adition to maternal and also will be added to paternal
#by default I would like to choose greed = 0
greed <- 0

#Explanation for the formula
#Basically, the formula use the modulus and quotient property of the alleles
#Before run the data using this script, alleles are converted in to several 'key values'
#Any key values can be used as long as the combination the key values raised for a new distinct value
#Quotient of modulus of my index value for the alleles give the weight for the maternal and paternal reads

for (i in 1:nallele){
  matsnp <- cbind(matsnp, snpdat$MAT%%(ploidy+1)^i)
}
for (i in 1:nallele){
  matsnp <- cbind(matsnp, matsnp[,i+1]%/%(ploidy+1)^(i-1))
}

for (i in 1:nallele){
  patsnp <- cbind(patsnp, snpdat$PAT%%(ploidy+1)^i)
}
for (i in 1:nallele){
  patsnp <- cbind(patsnp, patsnp[,i+1]%/%(ploidy+1)^(i-1))
}

#This sets all the number to 1 we don't want to have 2 in the calculation since it will multiplied the number of alleles
matsnp <- cbind(matsnp,ceiling(matsnp[,6:9]/(ploidy*2)))
patsnp <- cbind(patsnp,ceiling(patsnp[,6:9]/(ploidy*2)))

#greed variable take it role here (for explanation, see beginning of the script)
if (greed == 0){
	sumsnp <- abs(matsnp - patsnp)
} else {
	sumsnp <- matsnp + patsnp
}

newdmat <- sumsnp[,10:13]*matsnp[,10:13]
newdpat <- sumsnp[,10:13]*patsnp[,10:13]

newdmat <- cbind(newdmat, snpdat[,3:6]*newdmat)
newdpat <- cbind(newdpat, snpdat[,3:6]*newdpat)
newdmat <- cbind(newdmat, rowSums(newdmat[,1:4]))
names(newdmat)[9] <- "SUM"
newdpat <- cbind(newdpat, rowSums(newdpat[,1:4]))
names(newdpat)[9] <- "SUM"
newdmat <- cbind(newdmat, rowSums(newdmat[,5:8]))
names(newdmat)[10] <- "MATSUM"
newdpat <- cbind(newdpat, rowSums(newdpat[,5:8]))
names(newdpat)[10] <- "PATSUM"

sumtype <- rowSums(sumsnp[,10:13])
summat <- newdmat$MATSUM
sumpat <- newdpat$PATSUM
sumpar <- summat + sumpat

#Create another table
snpdat <- snpdat[,1:7]
snpdat <- cbind(snpdat, as.factor(genecat))
snpdat <- cbind(snpdat, as.factor(exoncat))
snpdat$POS[snpdat$POS > 0] <- 1
snpdat <- cbind(snpdat, as.numeric(sumtype))
snpdat <- cbind(snpdat, as.numeric(newdmat$MATSUM))
snpdat <- cbind(snpdat, as.numeric(newdpat$PATSUM))
snpdat <- cbind(snpdat, as.numeric(sumpar))
names(snpdat) <- c("NAMES","POS", "NUMA", "NUMC", "NUMG", "NUMT", "NUMN", "GENECAT", "EXON", "TYPE", "MATSUM", "PATSUM", "PARSUM")
snpdat <- snpdat[snpdat$POS > 0,]
snpdat <- snpdat[snpdat$TYPE > 1,]
snpdat <- snpdat[c("NAMES", "POS", "TYPE", "GENECAT", "EXON", "MATSUM","PATSUM","PARSUM")]

#This block is to identify maximum number of SNPs per GENECAT
snpsort <- snpdat[snpdat$GENECAT != "NA / NA / NA / NA",]
snpsort <- within(snpsort, RANK <- ave(-1 * PARSUM,GENECAT,FUN=rank, ties.method= "first"))
snpsort <- snpsort[order(snpsort$GENECAT),]
snpsort <- snpsort[c("GENECAT","RANK")]
snpcomp <- aggregate(.~ GENECAT, snpsort, max)

write.table(snpdat, file= opt$output, sep= "\t", quote= FALSE, row.names= FALSE, col.names= TRUE)
