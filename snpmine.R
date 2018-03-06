#!/usr/bin/env Rscript
##  MRH 26082013
##
##  snpmine.R v0, script to fill the SNP table with coverage
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
##	  
##  Parental file should be in tabular format with these columns: #Chromosome_number	Position_in_chromosome	Genome_Reference_Nucleotide	Genotype_in_N/N_format
##
##  input files:
##  2 filtered vcf file (one for maternal and the other one for paternal input)
##
##  output files:
##  compiled SNP list for both parents

print("For help: snpmine.R --help")

#The line exactly below is to clear all variables
rm(list=ls(all=TRUE))

#The line exactly below is to grab the arguments from the terminal command line
args <- commandArgs(TRUE)

#Need optparser and GenomicRanges Packages from Bioconductor

if(require('optparse')){
	print("optparse is available")
	} else {
	print("trying to install optparse")
	install.packages('optparse', repos = "http://R-Forge.R-project.org")
    if(require("optparse")){
		print("optparse is loaded")
    } else {
		stop("couldn't install optparse")
    }
}

library('optparse')

if(require('data.table')){
	print("data.table is available")
	} else {
	print("trying to install data.table")
	source("http://bioconductor.org/biocLite.R")
	biocLite('data.table', ask = FALSE)
	if(require("data.table")){
		print("data.table is loaded")
    } else {
		stop("couldn't install data.table")
    }
}

library('data.table')

option_list = list(
	make_option(c("-a", "--file1"), type = "character", default = NULL, help = "parent 1 path", metavar = "character"),
	make_option(c("-b", "--file2"), type = "character", default = NULL, help = "parent 2 path", metavar = "character"),
	make_option(c("-c", "--coverage1"), type = "character", default = NULL, help = "parent 1 coverage file", metavar = "character"),
	make_option(c("-d", "--coverage2"), type = "character", default = NULL, help = "parent 2 coverage file", metavar = "character"),
	make_option(c("-e", "--annotation"), type = "character", default = NULL, help = "dataset file name", metavar = "character"),
    make_option(c("-y", "--out1"), type = "character", default = NULL, help = "output file name", metavar = "character"),
    make_option(c("-z", "--out2"), type = "character", default = NULL, help = "output file name", metavar = "character"),
#this is a test feature
    make_option(c("-p", "--ploidy"), type="integer", default = 2, help = "input ploidy level (default = 2)", metavar = "integer")
)
 
opt_parser = OptionParser(option_list= option_list)
opt = parse_args(opt_parser)

if(require('GenomicRanges')){
	print("GenomicRanges is available")
	} else {
	print("trying to install GenomicRanges")
	source("http://bioconductor.org/biocLite.R")
	biocLite('GenomicRanges', ask = FALSE)
    if(require("GenomicRanges")){
        print("GenomicRanges is loaded")
    } else {
		stop("couldn't install GenomicRanges")
    }
}

library('GenomicRanges')

#Load SNP data. SNP table can been calculated manually using Libreoffice Calc or vcf-to-tab from vcftools
snp48 <- read.table(opt$file1, sep = '\t', header = FALSE)
snp75 <- read.table(opt$file2, sep = '\t', header = FALSE)
print("SNP data load success!")

names(snp48) <- c("CHROM","POS","REF","GENO")
snp48$SOURCE <- "cr48"
snp48 <- snp48[c("SOURCE", "CHROM", "POS", "REF", "GENO")]

names(snp75) <- c("CHROM","POS","REF","GENO")
snp75$SOURCE <- "cr75"
snp75 <- snp75[c("SOURCE", "CHROM", "POS", "REF", "GENO")]

#Split the 4th column (GENO) into two new column A1 and A2 and remove the first Column
#I don't know how to handle list created by lapply, so I need to convert them into character by using as.character
snp48$A1 <- lapply(strsplit(as.character(snp48$GENO), "/"),"[",1)
snp48$A2 <- lapply(strsplit(as.character(snp48$GENO), "/"),"[",2)
snp75$A1 <- lapply(strsplit(as.character(snp75$GENO), "/"),"[",1)
snp75$A2 <- lapply(strsplit(as.character(snp75$GENO), "/"),"[",2)
snp48$A1 <- as.character(snp48$A1)
snp48$A2 <- as.character(snp48$A2)
snp75$A1 <- as.character(snp75$A1)
snp75$A2 <- as.character(snp75$A2)

print(paste("Ploidy level: working with ", opt$ploidy, "n"))

##This block is an old block of deprected legacy function, soon it will be replaced with more sophisticated code

##Assign special number to specific allele
##evaluator will replace called genotype with specific number
evaluator <- function(x){
	if(x == "A"){return((opt$ploidy+1)^0)};
	if(x == "T"){return((opt$ploidy+1)^1)};
	if(x == "G"){return((opt$ploidy+1)^2)};
	if(x == "C"){return((opt$ploidy+1)^3)}
}
#evaluator2 is a function that in special case where the genotype was not called in other snp calling process will be replaced with the reference sequence genotype.
evaluator2 <- function(x){
	return(evaluator(x)*opt$ploidy)
}

#snp48$V1 <- lapply(snp48$A1, evaluator)
#snp48$V2 <- lapply(snp48$A2, evaluator)
#snp75$V1 <- lapply(snp75$A1, evaluator)
#snp75$V2 <- lapply(snp75$A2, evaluator)
#snp48$A1 <- NULL
#snp48$A2 <- NULL
#snp75$A1 <- NULL
#snp75$A2 <- NULL

#snp48$SMALLE <- as.numeric(snp48$V1) + as.numeric(snp48$V2)
#snp48$SMALLE <- as.character(snp48$SMALLE)
#snp75$SMALLE <- as.numeric(snp75$V1) + as.numeric(snp75$V2)
#snp75$SMALLE <- as.character(snp75$SMALLE)
#snp48$V1 <- NULL
#snp48$V2 <- NULL
#snp75$V1 <- NULL
#snp75$V2 <- NULL

##end of block

##This block below is an alternative way to convert allele to number, however, introduces bugs that haven't been tested.

snp48$SMALLE <- gsub("A", (2+1)^0, snp48$GENO)
snp48$SMALLE <- gsub("T", (2+1)^1, snp48$SMALLE)
snp48$SMALLE <- gsub("G", (2+1)^2, snp48$SMALLE)
snp48$SMALLE <- gsub("C", (2+1)^3, snp48$SMALLE)

snp75$SMALLE <- gsub("A", (2+1)^0, snp75$GENO)
snp75$SMALLE <- gsub("T", (2+1)^1, snp75$SMALLE)
snp75$SMALLE <- gsub("G", (2+1)^2, snp75$SMALLE)
snp75$SMALLE <- gsub("C", (2+1)^3, snp75$SMALLE)

calctable <- NULL
calculationtab <- NULL

calctable <- data.table(strsplit(as.character(snp48$SMALLE),'/',fixed=TRUE))
calctable$indexes <- seq(1,length(snp48$SMALLE),1)
calctable$num <- cbind(sapply(strsplit(as.character(calctable$V1),',',fixed=TRUE), length))
calculationtab <- data.frame(index= rep(calctable$indexes,calctable$num), value= as.numeric(unlist(sapply(calctable$V1, strsplit, ',', fixed=TRUE))))
calctable <- cbind(calctable, geno= aggregate(calculationtab$value, by= list(calculationtab$index), FUN= sum, na.rm = TRUE))

snp48$SMALLE <- calctable$geno.x

calctable <- NULL
calculationtab <- NULL

calctable <- data.table(strsplit(as.character(snp75$SMALLE),'/',fixed=TRUE))
calctable$indexes <- seq(1,length(snp75$SMALLE),1)
calctable$num <- cbind(sapply(strsplit(as.character(calctable$V1),',',fixed=TRUE), length))
calculationtab <- data.frame(index= rep(calctable$indexes,calctable$num), value= as.numeric(unlist(sapply(calctable$V1, strsplit, ',', fixed=TRUE))))
calctable <- cbind(calctable, geno= aggregate(calculationtab$value, by= list(calculationtab$index), FUN= sum, na.rm = TRUE))

snp75$SMALLE <- calctable$geno.x

calctable <- NULL
calculationtab <- NULL

##end of block

print("Indexing success!")

#Creating Ranges for subsetting
pos48 <- IRanges(snp48$POS, snp48$POS)
pos75 <- IRanges(snp75$POS, snp75$POS)
#Creating GenomicRanges for subsetting
gpos48 <- GRanges(snp48$CHROM, pos48)
gpos75 <- GRanges(snp75$CHROM, pos75)

#Read genome coverage
cov48 <- read.table(opt$coverage1, header = FALSE)
names(cov48) <- c("CHROM","POS","COV")
poscov48 <- IRanges(cov48$POS, cov48$POS)
gposcov48 <- GRanges(cov48$CHROM, poscov48)
cov75 <- read.table(opt$coverage2, header = FALSE)
names(cov75) <- c("CHROM","POS","COV")
poscov75 <- IRanges(cov75$POS, cov75$POS)
gposcov75 <- GRanges(cov75$CHROM, poscov75)

#Mining coverage
overlap48 <- findOverlaps(gpos48, gposcov48, select="first")
over <- as.matrix(overlap48)
#Subset the cov48 for position in snp48
cover48 <- cov48$COV[over]
newsnp48 <- cbind(snp48,cover48)
overlap48b <- findOverlaps(gpos48, gposcov75, select="first")
overb <- as.matrix(overlap48b)
#Subset the cov75 for position in snp48
cover75 <- cov75$COV[overb]
newsnp48 <- cbind(newsnp48,cover75)

#Zeroing NAs
newsnp48$cover75[is.na(newsnp48$cover75)] <- 0

#Mining coverage
overlap75 <- findOverlaps(gpos75, gposcov48,  select="first")
over <- as.matrix(overlap75)
#Subset the cov48 for position in snp75
cover48 <- cov48$COV[over]
newsnp75 <- cbind(snp75,cover48)
overlap75b <- findOverlaps(gpos75, gposcov75,  select="first")
overb <- as.matrix(overlap75b)
#Subset the cov75 for position in snp75
cover75 <- cov75$COV[overb]
newsnp75 <- cbind(newsnp75,cover75)

#Zeroing NAs
newsnp75$cover48[is.na(newsnp75$cover48)] <- 0

print("Coverage assignment success!")

#Find overlap between queries
overlapping <- findOverlaps(gpos48, gpos75)
overlapping <- as.matrix(overlapping)
newsnp48$ALT <- lapply(snp48$REF, evaluator2)
newsnp48$ALT <- as.character(newsnp48$ALT)
newsnp48$ALT[overlapping[,1]] <- snp75$SMALLE[overlapping[,2]]
newsnp48$OVER <- "FALSE"
newsnp48$OVER[overlapping[,1]] <- "TRUE"

overlapping <- findOverlaps(gpos75, gpos48)
overlapping <- as.matrix(overlapping)
newsnp75$ALT <- lapply(snp75$REF, evaluator2)
newsnp75$ALT <- as.character(newsnp75$ALT)
newsnp75$ALT[overlapping[,1]] <- snp48$SMALLE[overlapping[,2]]
newsnp75$OVER <- "FALSE"
newsnp75$OVER[overlapping[,1]] <- "TRUE"

#Read Annotation
annot <- read.table(opt$annotation, header = FALSE)
names(annot) <- c("seqid","source","type","start","end","score","strand","phase","attributes")
#Selecting exon annotation
exons <- annot[annot$type == 'exon',]
#Creating Ranges for subsetting
posexons <- IRanges(exons$start, exons$end)
#Creating GenomicRanges for subsetting
gposexons <- GRanges(exons$seqid, posexons, exons$strand, exons$type)

print("Annotation loaded!")

#Selecting gene annotation
genes <- annot[annot$type == 'gene',]
#Creating Ranges for subsetting
posgenes <- IRanges(genes$start, genes$end)
#Creating GenomicRanges for subsetting
gposgenes <- GRanges(genes$seqid, posgenes, genes$strand, genes$type)

#Find overlapping exons
overexons <- findOverlaps(gpos48, gposexons, select = "first")
over <- as.matrix(overexons)
exonhit <-  exons[over,]
exoncat <-  paste(exonhit$seqid,"/",exonhit$start,"/",exonhit$end,"/",exonhit$strand)
newsnp48 <- cbind(newsnp48,exoncat)

overexons <- findOverlaps(gpos75, gposexons, select = "first")
over <- as.matrix(overexons)
exonhit <-  exons[over,]
exoncat <-  paste(exonhit$seqid,"/",exonhit$start,"/",exonhit$end,"/",exonhit$strand)
newsnp75 <- cbind(newsnp75,exoncat)

print("Annotation assigment done!")

#Find overlapping genes
#DISTINCT field is a field which contain boolean information whether the parental genome alleles are the same or not. TRUE if they are different, FALSE if they are the same
overgenes <- findOverlaps(gpos48, gposgenes, select = "first")
over <- as.matrix(overgenes)
genehit <-  genes[over,]
genecat <-  paste(genehit$seqid,"/",genehit$start,"/",genehit$end,"/",genehit$strand)
newsnp48 <- cbind(newsnp48,genecat)
newsnp48$DISTINCT <- (newsnp48$SMALLE != newsnp48$ALT)

overgenes <- findOverlaps(gpos75, gposgenes, select = "first")
over <- as.matrix(overgenes)
genehit <-  genes[over,]
genecat <-  paste(genehit$seqid,"/",genehit$start,"/",genehit$end,"/",genehit$strand)
newsnp75 <- cbind(newsnp75,genecat)
newsnp75$DISTINCT <- (newsnp75$SMALLE != newsnp75$ALT)

#Filtering and Cosmetics

#filsnp48 <- newsnp48[newsnp48$cover48 <= 200 & newsnp48$cover75 >= 20 & newsnp48$cover75 <= 200,]
#Line below are added for anticipating unpredicted behaviour of stampy alignment. Remove (or put comments) to make it work
filsnp48 <- newsnp48[newsnp48$cover48 >= 20 & newsnp48$cover48 <= 200 & newsnp48$cover75 >= 20 & newsnp48$cover75 <= 200,]
filsnp48nov <- filsnp48[filsnp48$OVER == "FALSE",]
colnames(filsnp48)[which(names(filsnp48)=='SMALLE')] <- "MAT"
colnames(filsnp48)[which(names(filsnp48)=='ALT')] <- "PAT"

#filsnp75 <- newsnp75[newsnp75$cover75 <= 200 & newsnp75$cover48 >= 20 & newsnp75$cover48 <= 200,]
#Line below are added for anticipating unpredicted behaviour of stampy alignment. Remove (or put comments) to make it work
filsnp75 <- newsnp75[newsnp75$cover75 >= 20 & newsnp75$cover75 <= 200 & newsnp75$cover48 >= 20 & newsnp75$cover48 <= 200,]
filsnp75nov <- filsnp75[filsnp75$OVER == "FALSE",]
colnames(filsnp75)[which(names(filsnp75)=='SMALLE')] <- "MAT"
colnames(filsnp75)[which(names(filsnp75)=='ALT')] <- "PAT"

#Switching the allele to match to the maternal base
matsnpfilcomp <- rbind(filsnp48nov,filsnp75nov)
colnames(matsnpfilcomp)[which(names(matsnpfilcomp)=='SMALLE')] <- "MAT"
colnames(matsnpfilcomp)[which(names(matsnpfilcomp)=='ALT')] <- "PAT"
#SMALLE column in snp75 is its alllele which is paternal from Cr48 x Cr75 side and otherwise for ALT
matsnpfilcomp$MAT[matsnpfilcomp$SOURCE == "cr75"] <- filsnp75nov$ALT
matsnpfilcomp$PAT[matsnpfilcomp$SOURCE == "cr75"] <- filsnp75nov$SMALLE

#Switching the allele to match to the maternal base
patsnpfilcomp <- rbind(filsnp48nov,filsnp75nov)
colnames(patsnpfilcomp)[which(names(patsnpfilcomp)=='SMALLE')] <- "MAT"
colnames(patsnpfilcomp)[which(names(patsnpfilcomp)=='ALT')] <- "PAT"
#SMALLE column in snp75 is its alllele which is paternal from Cr75 x Cr48 side and otherwise for ALT
patsnpfilcomp$MAT[patsnpfilcomp$SOURCE == "cr48"] <- filsnp48nov$ALT
patsnpfilcomp$PAT[patsnpfilcomp$SOURCE == "cr48"] <- filsnp48nov$SMALLE

print("Switching succes!")

#If Distinct is TRUE then the sequence must be overlapped and has different value which is needed for analysis. However, overlapping region with same value automatically discarded
filsnp48 <- filsnp48[filsnp48$OVER == "TRUE",]
filsnp48 <- filsnp48[filsnp48$DISTINCT == "TRUE",]
matsnpfilcomp <- rbind(matsnpfilcomp, filsnp48)
#add/remove (comment) line below to make the downstream analysis faster, however, this process might reduce certain possibilities in assessing region without annotation
matsnpfilcomp <- matsnpfilcomp[matsnpfilcomp$genecat != "NA / NA / NA / NA",]


filsnp75 <- filsnp75[filsnp75$OVER == "TRUE",]
filsnp75 <- filsnp75[filsnp75$DISTINCT == "TRUE",]
patsnpfilcomp <- rbind(patsnpfilcomp, filsnp75)
#add/remove (comment) line below to make the downstream analysis faster, however, this process might reduce certain possibilities in assessing region without annotation
patsnpfilcomp <- patsnpfilcomp[patsnpfilcomp$genecat != "NA / NA / NA / NA",]

write.table(matsnpfilcomp,opt$out1, row.names=FALSE, sep="\t", quote=FALSE)
write.table(patsnpfilcomp,opt$out2, row.names=FALSE, sep="\t", quote=FALSE)

print("Overlap analysis done!")
print("All analysis is finished!")
