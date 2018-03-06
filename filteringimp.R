#MRH 18112013
#!/usr/bin/env Rscript
#MRH 18112013
##  
##  filteringimp.R v0, This script is used to filter imprinting genes based on their consistency in two replicates.
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


#Opt section
library('optparse')
option_list = list(
	make_option(c("-a", "--input1"), type = "character", default = NULL, help = "parent 1", metavar = "character"),
	make_option(c("-b", "--input2"), type = "character", default = NULL, help = "parent 2", metavar = "character"),
	make_option(c("-c", "--compile"), type = "character", default = NULL, help = "compiled output file", metavar = "character"),
	make_option(c("-r", "--rawoutput"), type = "character", default = NULL, help = "raw output file", metavar = "character"),
	make_option(c("-o", "--output"), type = "character", default = NULL, help = "compiled output file", metavar = "character")
)

opt_parser = OptionParser(option_list= option_list)
opt = parse_args(opt_parser)
#End of opt section

rep1 <- read.table(opt$input1, header= TRUE, sep= "\t")
#rep1 <- rep1[c("NAMES", "POS", "MAT", "PAT", "summat", "sumpat", "parchitest", "mattwobintest", "pattwobintest")]
names(rep1) = paste(names(rep1),sep='',1)
rep2 <- read.table(opt$input2, header= TRUE, sep= "\t")
names(rep2) = paste(names(rep2),sep='',2)
#rep2 <- rep2[c("NAMES", "POS", "MAT", "PAT", "summat", "sumpat", "parchitest", "mattwobintest", "pattwobintest")]
contain <- merge (rep1, rep2, by.x= "NAMES1", by.y= "NAMES2", all= TRUE)
colnames(contain)[1] <- "NAMES"

#NA/NA/NA/NA Removal Exon and Genecat
binder <- data.frame(as.character(contain$EXON1),as.character(contain$EXON2), stringsAsFactors= FALSE)
names(binder) <- c('EXON1','EXON2')
binder$EXON1[is.na(binder$EXON1)] <- 'NA / NA / NA / NA'
binder$EXON2[is.na(binder$EXON2)] <- 'NA / NA / NA / NA'
list1 <- as.vector(binder$EXON1)
list2 <- as.vector(binder$EXON2)
binder$EXON1 <- ifelse(binder$EXON1 == 'NA / NA / NA / NA', list2, list1)
binder$EXON2 <- ifelse(binder$EXON2 == 'NA / NA / NA / NA', list1, list2)
contain$EXON1 <- binder$EXON1
contain$EXON2 <- binder$EXON2

binder <- data.frame(as.character(contain$GENECAT1),as.character(contain$GENECAT2), stringsAsFactors= FALSE)
names(binder) <- c('GENECAT1','GENECAT2')
binder$GENECAT1[is.na(binder$GENECAT1)] <- 'NA / NA / NA / NA'
binder$GENECAT2[is.na(binder$GENECAT2)] <- 'NA / NA / NA / NA'
list1 <- as.vector(binder$GENECAT1)
list2 <- as.vector(binder$GENECAT2)
binder$GENECAT1 <- ifelse(binder$GENECAT1 == 'NA / NA / NA / NA', list2, list1)
binder$GENECAT2 <- ifelse(binder$GENECAT2 == 'NA / NA / NA / NA', list1, list2)
contain$GENECAT1 <- binder$GENECAT1
contain$GENECAT2 <- binder$GENECAT2

#The CORRECTIONBLOCK. The calculation of MRATIO has been done in *det.csv file. However, some information might be gone. While the calculation below seems redundant, it is needed to recover missing value from merging table.
#contain$sumpar1 <- contain$sumpat1 + contain$summat1
#contain$sumpar2 <- contain$sumpat2 + contain$summat2

#contain$MRATIO1 <- contain$MATSUM1/(2*contain$PATSUM1 + contain$MATSUM1)
#contain$MRATIO2 <- contain$MATSUM2/(2*contain$PATSUM2 + contain$MATSUM2)
elselist1 <- as.vector(contain$MATSUM1/(2*contain$PATSUM1 + contain$MATSUM1))
elselist2 <- as.vector(contain$MRATIO1)
contain$MRATIO1 <- ifelse(is.na(contain$MRATIO1), elselist1, elselist2)

elselist1 <- as.vector(contain$MATSUM2/(2*contain$PATSUM2 + contain$MATSUM2))
elselist2 <- as.vector(contain$MRATIO2)
contain$MRATIO2 <- ifelse(is.na(contain$MRATIO2), elselist1, elselist2)
#End of CORRECTIONBLOCK

contain$IMP1[contain$MRATIO1 == 0.5] <- "N"
contain$IMP1[contain$MRATIO1 > 0.5] <- "M"
contain$IMP1[contain$MRATIO1 < 0.5] <- "P"
contain$IMP2[contain$MRATIO2 == 0.5] <- "N"
contain$IMP2[contain$MRATIO2 > 0.5] <- "M"
contain$IMP2[contain$MRATIO2 < 0.5] <- "P"

contain$CONS <- contain$IMP1 == contain$IMP2
contain$PARFRAC1 <- apply(cbind(1-contain$MRATIO1, contain$MRATIO1),1,max)
contain$PARFRAC2 <- apply(cbind(1-contain$MRATIO2, contain$MRATIO2),1,max)

#Evidence for non-imprinted genes
notimprint <- subset(contain, PARFRAC1<= 0.8 & PARFRAC2<= 0.8)
notimprint <- subset(notimprint, TOTSUM1>= 10 & TOTSUM2>= 10)
notimprint <- subset(notimprint, genchitest1>= 0.05 & genchitest2>= 0.05)
#write.table(notimprint, file= "/home/rocky/Imprinting/notimprint2.csv", sep= "\t", quote= FALSE, row.names= FALSE, col.names= TRUE)

#Detection of Dominance Effect
contain$DOMINANCE <- ifelse(contain$CONS == FALSE & (contain$genchitest1 <= 0.05 | contain$genchitest2 <= 0.05), TRUE, FALSE)

#Detection of Parental Bias
contain$BIAS <- ifelse(contain$CONS == TRUE & ((contain$genchitest1 <= 0.05 & contain$genchitest2 > 0.05) | (contain$genchitest2 <= 0.05 & contain$genchitest1 > 0.05)), TRUE, FALSE)
parbias <- contain[contain$BIAS == TRUE,]
parbias <- subset(parbias, PARFRAC1>= 0.8 | PARFRAC2>= 0.8)

#Conservative way
#parbias <- subset(parbias, PARSUM1>= 20 & PARSUM2>= 20)
#Greedy way
parbias <- subset(parbias, TOTSUM1>= 10 & TOTSUM2>= 10)

nrow(parbias[parbias$IMP1 == "M",])
nrow(parbias[parbias$IMP1 == "P",])
parbias$IMPTYPE <- "B"
#write.table(parbias, file= "/home/rocky/Imprinting/filparbiasr2.csv", sep= "\t", quote= FALSE, row.names= FALSE, col.names= TRUE)

#Filtering p-val and consistancy
report <- subset(contain, CONS == TRUE)

#Filter PARFRAC if over than 4 fold than normal expression
report <- subset(report, PARFRAC1>= 0.8 & PARFRAC2>= 0.8)

#Conservative way
#report <- subset(report, PARSUM1>= 20 & PARSUM2>= 20)
#Greedy way
report <- subset(report, TOTSUM1>= 10 & TOTSUM2>= 10)
report <- subset(report, genchitest1<= 0.01 & genchitest2<= 0.01)

nrow(report[report$IMP1 == "M",])
nrow(report[report$IMP1 == "P",])
weakimp <- report

#Filter PARFRAC if over than 8 fold than normal expression
report <- subset(report, PARFRAC1>= 0.888889 & PARFRAC2>= 0.888889)
nrow(report[report$IMP1 == "M",])
nrow(report[report$IMP1 == "P",])
moderateimp <- report

#Filter PARFRAC if over than 16 fold than normal expression
report <- subset(report, PARFRAC1>= 0.941176 & PARFRAC2>= 0.941176)
nrow(report[report$IMP1 == "M",])
nrow(report[report$IMP1 == "P",])
strongimp <- report

weakimp <- weakimp[!weakimp$NAMES %in% moderateimp$NAMES,]
weakimp$IMPTYPE <- "W"
#write.table(weakimp, file= "/home/rocky/Imprinting/filimprintwear2.csv", sep= "\t", quote= FALSE, row.names= FALSE, col.names= TRUE)
moderateimp <- moderateimp[!moderateimp$NAMES %in% strongimp$NAMES,]
moderateimp$IMPTYPE <- "M"
#write.table(moderateimp, file= "/home/rocky/Imprinting/filimprintmodr2.csv", sep= "\t", quote= FALSE, row.names= FALSE, col.names= TRUE)
strongimp$IMPTYPE <- "S"

newtable <- rbind(parbias, weakimp, moderateimp, strongimp)
#write.table(strongimp, file= "/home/rocky/Imprinting/filimprintstrongr2.csv", sep= "\t", quote= FALSE, row.names= FALSE, col.names= TRUE)

write.table(contain, file= opt$rawoutput, sep= "\t", quote= FALSE, row.names= FALSE, col.names= TRUE)
write.table(newtable, file= opt$compile, sep= "\t", quote= FALSE, row.names= FALSE, col.names= TRUE)

aliastab <- data.frame(GENECAT= unique(c(as.vector(contain$GENECAT1),as.vector(contain$GENECAT2))))
subtab <- data.frame(GENECAT= contain$GENECAT1, MRATIO1= contain$MRATIO1, TOTSUM1= contain$TOTSUM1, genchitest1= contain$genchitest1, MRATIO2= contain$MRATIO2,  TOTSUM2= contain$TOTSUM2, genchitest2= contain$genchitest2)
newtab <- merge(aliastab, aggregate(MRATIO1 ~ GENECAT, subtab, FUN= max, na.rm=TRUE), all.x= TRUE)
newtab <- merge(newtab, aggregate(TOTSUM1 ~ GENECAT, subtab, FUN= max, na.rm=TRUE), all.x= TRUE)
newtab <- merge(newtab, aggregate(genchitest1 ~ GENECAT, subtab, FUN= max, na.rm=TRUE), all.x= TRUE)
newtab <- merge(newtab, aggregate(MRATIO2 ~ GENECAT, subtab, FUN= max, na.rm=TRUE), all.x= TRUE)
newtab <- merge(newtab, aggregate(TOTSUM2 ~ GENECAT, subtab, FUN= max, na.rm=TRUE), all.x= TRUE)
newtab <- merge(newtab, aggregate(genchitest2 ~ GENECAT, subtab, FUN= max, na.rm=TRUE), all.x= TRUE)

write.table(newtab, file= paste("simplified",opt$rawoutput, sep= ""), sep= "\t", quote= FALSE, row.names= FALSE, col.names= TRUE)

#To help in filtering and compiling
compiling <- newtable[c("GENECAT1","IMP1", "IMPTYPE")]
nrow(compiling)
compiling <- compiling[!duplicated(compiling),]
names(compiling) <- c("GENECAT","IMP", "IMPTYPE")
nrow(compiling)
mcompiling <- compiling[compiling$IMP == 'M',]
nrow(mcompiling[mcompiling$IMPTYPE == 'B',])
nrow(mcompiling[mcompiling$IMPTYPE == 'W',])
nrow(mcompiling[mcompiling$IMPTYPE == 'M',])
nrow(mcompiling[mcompiling$IMPTYPE == 'S',])
pcompiling <- compiling[compiling$IMP == 'P',]
nrow(pcompiling[pcompiling$IMPTYPE == 'B',])
nrow(pcompiling[pcompiling$IMPTYPE == 'W',])
nrow(pcompiling[pcompiling$IMPTYPE == 'M',])
nrow(pcompiling[pcompiling$IMPTYPE == 'S',])

filtertab <- read.table("filter.csv", header= TRUE, sep= "\t")
mcompiling <- mcompiling[mcompiling$IMPTYPE != 'B',]
mcompiling <- mcompiling[mcompiling$GENECAT %in% filtertab[,1],]
pcompiling <- pcompiling[pcompiling$IMPTYPE != 'B',]
filcompiling <- rbind(mcompiling, pcompiling)

write.table(compiling, file= opt$output, sep= "\t", quote= FALSE, row.names= FALSE, col.names= TRUE)
write.table(filcompiling, file= paste("filcom",opt$output, sep= ""), sep= "\t", quote= FALSE, row.names= FALSE, col.names= TRUE)

