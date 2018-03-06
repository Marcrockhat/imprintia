#!/usr/bin/env Rscript
#MRH 18112013
##  
##  stattest.R v0, This script is used to identify maternal and paternal SNP and calculate the p-value on proportion test based on several methods.
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

#library('WRS')

args <- commandArgs(TRUE)

#Opt section
library('optparse')
option_list = list(
	make_option(c("-i", "--input"), type = "character", default = NULL, help = "parent input", metavar = "character"),
	make_option(c("-o", "--output"), type = "character", default = NULL, help = "output file", metavar = "character")
)

opt_parser = OptionParser(option_list= option_list)
opt = parse_args(opt_parser)
#End of opt section

snpdat <- read.table(opt$input, header= TRUE, sep= "\t")

#Filter for not useable reads. The reason for this is: it was estimated that one genes contain 10 SNPs. It can raise bias for low covered sequences where the SNP might come from a couple number of reads.
#Depends on situation, turn this on and off
#snpdat <- snpdat[snpdat$PARSUM >= 20,]

parchitest <- NULL

summat <- snpdat$MATSUM
sumpat <- snpdat$PATSUM
sumpar <- snpdat$PARSUM

#Test the chisquare test
for (i in 1:length(sumpar)){
  if(summat[i] > 0 | sumpat[i] > 0){
    chitest <- chisq.test(c(summat[i],sumpat[i]), p= c(2/3, 1/3))
    parchitest <- rbind(parchitest, chitest$p.value)
  } else {
    parchitest <- rbind(parchitest, NA)
  }
}

locdat <- snpdat
locdat <- cbind(locdat, as.vector(p.adjust(parchitest, "BH")))
names(locdat)[9] <- 'locchitest'

exdat <- locdat[locdat$EXON != "NA / NA / NA / NA",]
exdat <- aggregate(. ~ EXON, data= exdat, FUN= sum)

summat <- exdat$MATSUM
sumpat <- exdat$PATSUM
sumpar <- exdat$PARSUM

parchitest <- NULL

#Second Block for Per Exon Test
#Test the maternal binomial proportion test

#Test the chisquare test
for (i in 1:length(sumpar)){
  if(summat[i] > 0 | sumpat[i] > 0){
    chitest <- chisq.test(c(summat[i],sumpat[i]), p= c(2/3, 1/3))
    parchitest <- rbind(parchitest, chitest$p.value)
  } else {
    parchitest <- rbind(parchitest, NA)
  }
}

exdat <- exdat['EXON']
exdat <- cbind(exdat, as.vector(p.adjust(parchitest, "BH")))
names(exdat)[2] <- 'exchitest'
exdatmer <- merge(locdat, exdat, by.x= 'EXON', by.y= 'EXON', all.x= TRUE)

gendat <- snpdat[snpdat$GENECAT != "NA / NA / NA / NA",]
gendat <- aggregate(. ~ GENECAT, data= gendat, FUN= sum)

summat <- gendat$MATSUM
sumpat <- gendat$PATSUM
sumpar <- gendat$PARSUM

parchitest <- NULL

#Test the chisquare test
for (i in 1:length(sumpar)){
  if(summat[i] > 0 | sumpat[i] > 0){
    chitest <- chisq.test(c(summat[i],sumpat[i]), p= c(2/3, 1/3))
    parchitest <- rbind(parchitest, chitest$p.value)
  } else {
    parchitest <- rbind(parchitest, NA)
  }
}

gendat <- gendat['GENECAT']
gendat <- cbind(gendat, as.vector(p.adjust(parchitest, "BH")))
names(gendat)[2] <- 'genchitest'

#This line specially applied to gendat
gendat$MRATIO <- summat/(2*sumpat + summat)
gendat$TOTSUM <- sumpar

tempdat <- merge(exdatmer, gendat, by.x= 'GENECAT', by.y= 'GENECAT', all.x= TRUE)
tempdat <- tempdat[c('NAMES', 'EXON', 'GENECAT',  'POS', 'TYPE', 'MATSUM', 'PATSUM', 'PARSUM', 'locchitest', 'exchitest', 'genchitest', 'MRATIO', 'TOTSUM')]
tempdat <- tempdat[order(tempdat$NAMES),]

write.table(tempdat, file= opt$output, sep= "\t",quote= FALSE, row.names= FALSE, col.names= TRUE)
