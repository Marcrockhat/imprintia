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

#If you are running imprintia as a pipeline using imprintia.sh. You will need to change the default value manually. Automatization will be added later to accomodate easy deployment.
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
snp_par_A <- read.table(opt$file1, sep = '\t', header = FALSE)
snp_par_B <- read.table(opt$file2, sep = '\t', header = FALSE)
print("SNP data load success!")

names(snp_par_A) <- c("CHROM","POS","REF","GENO")
snp_par_A$SOURCE <- "Par_A_Geno"
snp_par_A <- snp_par_A[c("SOURCE", "CHROM", "POS", "REF", "GENO")]

names(snp_par_B) <- c("CHROM","POS","REF","GENO")
snp_par_B$SOURCE <- "Par_B_Geno"
snp_par_B <- snp_par_B[c("SOURCE", "CHROM", "POS", "REF", "GENO")]

#Split the 4th column (GENO) into two new column A1 and A2 and remove the first Column
#I don't know how to handle list created by lapply, so I need to convert them into character by using as.character
snp_par_A$A1 <- lapply(strsplit(as.character(snp_par_A$GENO), "/"),"[",1)
snp_par_A$A2 <- lapply(strsplit(as.character(snp_par_A$GENO), "/"),"[",2)
snp_par_B$A1 <- lapply(strsplit(as.character(snp_par_B$GENO), "/"),"[",1)
snp_par_B$A2 <- lapply(strsplit(as.character(snp_par_B$GENO), "/"),"[",2)
snp_par_A$A1 <- as.character(snp_par_A$A1)
snp_par_A$A2 <- as.character(snp_par_A$A2)
snp_par_B$A1 <- as.character(snp_par_B$A1)
snp_par_B$A2 <- as.character(snp_par_B$A2)

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

#snp_par_A$V1 <- lapply(snp_par_A$A1, evaluator)
#snp_par_A$V2 <- lapply(snp_par_A$A2, evaluator)
#snp_par_B$V1 <- lapply(snp_par_B$A1, evaluator)
#snp_par_B$V2 <- lapply(snp_par_B$A2, evaluator)
#snp_par_A$A1 <- NULL
#snp_par_A$A2 <- NULL
#snp_par_B$A1 <- NULL
#snp_par_B$A2 <- NULL

#snp_par_A$SMALLE <- as.numeric(snp_par_A$V1) + as.numeric(snp_par_A$V2)
#snp_par_A$SMALLE <- as.character(snp_par_A$SMALLE)
#snp_par_B$SMALLE <- as.numeric(snp_par_B$V1) + as.numeric(snp_par_B$V2)
#snp_par_B$SMALLE <- as.character(snp_par_B$SMALLE)
#snp_par_A$V1 <- NULL
#snp_par_A$V2 <- NULL
#snp_par_B$V1 <- NULL
#snp_par_B$V2 <- NULL

##end of block

##This block below is an alternative way to convert allele to number, however, introduces bugs that haven't been tested.

snp_par_A$SMALLE <- gsub("A", (2+1)^0, snp_par_A$GENO)
snp_par_A$SMALLE <- gsub("T", (2+1)^1, snp_par_A$SMALLE)
snp_par_A$SMALLE <- gsub("G", (2+1)^2, snp_par_A$SMALLE)
snp_par_A$SMALLE <- gsub("C", (2+1)^3, snp_par_A$SMALLE)

snp_par_B$SMALLE <- gsub("A", (2+1)^0, snp_par_B$GENO)
snp_par_B$SMALLE <- gsub("T", (2+1)^1, snp_par_B$SMALLE)
snp_par_B$SMALLE <- gsub("G", (2+1)^2, snp_par_B$SMALLE)
snp_par_B$SMALLE <- gsub("C", (2+1)^3, snp_par_B$SMALLE)

calctable <- NULL
calculationtab <- NULL

calctable <- data.table(strsplit(as.character(snp_par_A$SMALLE),'/',fixed=TRUE))
calctable$indexes <- seq(1,length(snp_par_A$SMALLE),1)
calctable$num <- cbind(sapply(strsplit(as.character(calctable$V1),',',fixed=TRUE), length))
calculationtab <- data.frame(index= rep(calctable$indexes,calctable$num), value= as.numeric(unlist(sapply(calctable$V1, strsplit, ',', fixed=TRUE))))
calctable <- cbind(calctable, geno= aggregate(calculationtab$value, by= list(calculationtab$index), FUN= sum, na.rm = TRUE))

snp_par_A$SMALLE <- calctable$geno.x

calctable <- NULL
calculationtab <- NULL

calctable <- data.table(strsplit(as.character(snp_par_B$SMALLE),'/',fixed=TRUE))
calctable$indexes <- seq(1,length(snp_par_B$SMALLE),1)
calctable$num <- cbind(sapply(strsplit(as.character(calctable$V1),',',fixed=TRUE), length))
calculationtab <- data.frame(index= rep(calctable$indexes,calctable$num), value= as.numeric(unlist(sapply(calctable$V1, strsplit, ',', fixed=TRUE))))
calctable <- cbind(calctable, geno= aggregate(calculationtab$value, by= list(calculationtab$index), FUN= sum, na.rm = TRUE))

snp_par_B$SMALLE <- calctable$geno.x

calctable <- NULL
calculationtab <- NULL

##end of block

print("Indexing success!")

#Creating Ranges for subsetting
pos_par_A <- IRanges(snp_par_A$POS, snp_par_A$POS)
pos_par_B <- IRanges(snp_par_B$POS, snp_par_B$POS)
#Creating GenomicRanges for subsetting
gpos_par_A <- GRanges(snp_par_A$CHROM, pos_par_A)
gpos_par_B <- GRanges(snp_par_B$CHROM, pos_par_B)

#Read genome coverage
cov_par_A <- read.table(opt$coverage1, header = FALSE)
names(cov_par_A) <- c("CHROM","POS","COV")
poscov_par_A <- IRanges(cov_par_A$POS, cov_par_A$POS)
gposcov_par_A <- GRanges(cov_par_A$CHROM, poscov_par_A)
cov_par_B <- read.table(opt$coverage2, header = FALSE)
names(cov_par_B) <- c("CHROM","POS","COV")
poscov_par_B <- IRanges(cov_par_B$POS, cov_par_B$POS)
gposcov_par_B <- GRanges(cov_par_B$CHROM, poscov_par_B)

#Mining coverage
overlap_par_A <- findOverlaps(gpos_par_A, gposcov_par_A, select="first")
over <- as.matrix(overlap_par_A)
#Subset the cov_par_A for position in snp_par_A
cover_par_A <- cov_par_A$COV[over]
newsnp_par_A <- cbind(snp_par_A,cover_par_A)
overlap_par_Ab <- findOverlaps(gpos_par_A, gposcov_par_B, select="first")
overb <- as.matrix(overlap_par_Ab)
#Subset the cov_par_B for position in snp_par_A
cover_par_B <- cov_par_B$COV[overb]
newsnp_par_A <- cbind(newsnp_par_A,cover_par_B)

#Zeroing NAs
newsnp_par_A$cover_par_B[is.na(newsnp_par_A$cover_par_B)] <- 0

#Mining coverage
overlap_par_B <- findOverlaps(gpos_par_B, gposcov_par_A,  select="first")
over <- as.matrix(overlap_par_B)
#Subset the cov_par_A for position in snp_par_B
cover_par_A <- cov_par_A$COV[over]
newsnp_par_B <- cbind(snp_par_B,cover_par_A)
overlap_par_Bb <- findOverlaps(gpos_par_B, gposcov_par_B,  select="first")
overb <- as.matrix(overlap_par_Bb)
#Subset the cov_par_B for position in snp_par_B
cover_par_B <- cov_par_B$COV[overb]
newsnp_par_B <- cbind(newsnp_par_B,cover_par_B)

#Zeroing NAs
newsnp_par_B$cover_par_A[is.na(newsnp_par_B$cover_par_A)] <- 0

print("Coverage assignment success!")

#Find overlap between queries
overlapping <- findOverlaps(gpos_par_A, gpos_par_B)
overlapping <- as.matrix(overlapping)
newsnp_par_A$ALT <- lapply(snp_par_A$REF, evaluator2)
newsnp_par_A$ALT <- as.character(newsnp_par_A$ALT)
newsnp_par_A$ALT[overlapping[,1]] <- snp_par_B$SMALLE[overlapping[,2]]
newsnp_par_A$OVER <- "FALSE"
newsnp_par_A$OVER[overlapping[,1]] <- "TRUE"

overlapping <- findOverlaps(gpos_par_B, gpos_par_A)
overlapping <- as.matrix(overlapping)
newsnp_par_B$ALT <- lapply(snp_par_B$REF, evaluator2)
newsnp_par_B$ALT <- as.character(newsnp_par_B$ALT)
newsnp_par_B$ALT[overlapping[,1]] <- snp_par_A$SMALLE[overlapping[,2]]
newsnp_par_B$OVER <- "FALSE"
newsnp_par_B$OVER[overlapping[,1]] <- "TRUE"

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
overexons <- findOverlaps(gpos_par_A, gposexons, select = "first")
over <- as.matrix(overexons)
exonhit <-  exons[over,]
exoncat <-  paste(exonhit$seqid,"/",exonhit$start,"/",exonhit$end,"/",exonhit$strand)
newsnp_par_A <- cbind(newsnp_par_A,exoncat)

overexons <- findOverlaps(gpos_par_B, gposexons, select = "first")
over <- as.matrix(overexons)
exonhit <-  exons[over,]
exoncat <-  paste(exonhit$seqid,"/",exonhit$start,"/",exonhit$end,"/",exonhit$strand)
newsnp_par_B <- cbind(newsnp_par_B,exoncat)

print("Annotation assigment done!")

#Find overlapping genes
#DISTINCT field is a field which contain boolean information whether the parental genome alleles are the same or not. TRUE if they are different, FALSE if they are the same
overgenes <- findOverlaps(gpos_par_A, gposgenes, select = "first")
over <- as.matrix(overgenes)
genehit <-  genes[over,]
genecat <-  paste(genehit$seqid,"/",genehit$start,"/",genehit$end,"/",genehit$strand)
newsnp_par_A <- cbind(newsnp_par_A,genecat)
newsnp_par_A$DISTINCT <- (newsnp_par_A$SMALLE != newsnp_par_A$ALT)

overgenes <- findOverlaps(gpos_par_B, gposgenes, select = "first")
over <- as.matrix(overgenes)
genehit <-  genes[over,]
genecat <-  paste(genehit$seqid,"/",genehit$start,"/",genehit$end,"/",genehit$strand)
newsnp_par_B <- cbind(newsnp_par_B,genecat)
newsnp_par_B$DISTINCT <- (newsnp_par_B$SMALLE != newsnp_par_B$ALT)

#Filtering and Cosmetics

#filsnp_par_A <- newsnp_par_A[newsnp_par_A$cover_par_A <= 200 & newsnp_par_A$cover_par_B >= 20 & newsnp_par_A$cover_par_B <= 200,]
#Line below are added for anticipating unpredicted behaviour of stampy alignment. Remove (or put comments) to make it work
filsnp_par_A <- newsnp_par_A[newsnp_par_A$cover_par_A >= 20 & newsnp_par_A$cover_par_A <= 200 & newsnp_par_A$cover_par_B >= 20 & newsnp_par_A$cover_par_B <= 200,]
filsnp_par_Anov <- filsnp_par_A[filsnp_par_A$OVER == "FALSE",]
colnames(filsnp_par_A)[which(names(filsnp_par_A)=='SMALLE')] <- "MAT"
colnames(filsnp_par_A)[which(names(filsnp_par_A)=='ALT')] <- "PAT"

#filsnp_par_B <- newsnp_par_B[newsnp_par_B$cover_par_B <= 200 & newsnp_par_B$cover_par_A >= 20 & newsnp_par_B$cover_par_A <= 200,]
#Line below are added for anticipating unpredicted behaviour of stampy alignment. Remove (or put comments) to make it work
filsnp_par_B <- newsnp_par_B[newsnp_par_B$cover_par_B >= 20 & newsnp_par_B$cover_par_B <= 200 & newsnp_par_B$cover_par_A >= 20 & newsnp_par_B$cover_par_A <= 200,]
filsnp_par_Bnov <- filsnp_par_B[filsnp_par_B$OVER == "FALSE",]
colnames(filsnp_par_B)[which(names(filsnp_par_B)=='SMALLE')] <- "MAT"
colnames(filsnp_par_B)[which(names(filsnp_par_B)=='ALT')] <- "PAT"

#Switching the allele to match to the maternal base
matsnpfilcomp <- rbind(filsnp_par_Anov,filsnp_par_Bnov)
colnames(matsnpfilcomp)[which(names(matsnpfilcomp)=='SMALLE')] <- "MAT"
colnames(matsnpfilcomp)[which(names(matsnpfilcomp)=='ALT')] <- "PAT"
#SMALLE column in snp_par_B is its alllele which is paternal from Par_A_Geno x Par_B_Geno side and otherwise for ALT
matsnpfilcomp$MAT[matsnpfilcomp$SOURCE == "Par_B_Geno"] <- filsnp_par_Bnov$ALT
matsnpfilcomp$PAT[matsnpfilcomp$SOURCE == "Par_B_Geno"] <- filsnp_par_Bnov$SMALLE

#Switching the allele to match to the maternal base
patsnpfilcomp <- rbind(filsnp_par_Anov,filsnp_par_Bnov)
colnames(patsnpfilcomp)[which(names(patsnpfilcomp)=='SMALLE')] <- "MAT"
colnames(patsnpfilcomp)[which(names(patsnpfilcomp)=='ALT')] <- "PAT"
#SMALLE column in snp_par_B is its alllele which is paternal from Par_B_Geno x Par_A_Geno side and otherwise for ALT
patsnpfilcomp$MAT[patsnpfilcomp$SOURCE == "Par_A_Geno"] <- filsnp_par_Anov$ALT
patsnpfilcomp$PAT[patsnpfilcomp$SOURCE == "Par_A_Geno"] <- filsnp_par_Anov$SMALLE

print("Switching succes!")

#If Distinct is TRUE then the sequence must be overlapped and has different value which is needed for analysis. However, overlapping region with same value automatically discarded
filsnp_par_A <- filsnp_par_A[filsnp_par_A$OVER == "TRUE",]
filsnp_par_A <- filsnp_par_A[filsnp_par_A$DISTINCT == "TRUE",]
matsnpfilcomp <- rbind(matsnpfilcomp, filsnp_par_A)
#add/remove (comment) line below to make the downstream analysis faster, however, this process might reduce certain possibilities in assessing region without annotation
matsnpfilcomp <- matsnpfilcomp[matsnpfilcomp$genecat != "NA / NA / NA / NA",]


filsnp_par_B <- filsnp_par_B[filsnp_par_B$OVER == "TRUE",]
filsnp_par_B <- filsnp_par_B[filsnp_par_B$DISTINCT == "TRUE",]
patsnpfilcomp <- rbind(patsnpfilcomp, filsnp_par_B)
#add/remove (comment) line below to make the downstream analysis faster, however, this process might reduce certain possibilities in assessing region without annotation
patsnpfilcomp <- patsnpfilcomp[patsnpfilcomp$genecat != "NA / NA / NA / NA",]

write.table(matsnpfilcomp,opt$out1, row.names=FALSE, sep="\t", quote=FALSE)
write.table(patsnpfilcomp,opt$out2, row.names=FALSE, sep="\t", quote=FALSE)

print("Overlap analysis done!")
print("All analysis is finished!")
