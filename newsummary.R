#MRH 24052014
##
##  newsummary.R v0, compiled SNP list for both parents#Script to summarize imprintin and filter out not specific endospermn expressed genes
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
	make_option(c("-a", "--input1"), type = "character", default = NULL, help = "first replicate path", metavar = "character"),
	make_option(c("-b", "--input2"), type = "character", default = NULL, help = "second replicate path", metavar = "character"),
	make_option(c("-c", "--raw1"), type = "character", default = NULL, help = "first replicate path", metavar = "character"),
	make_option(c("-d", "--raw2"), type = "character", default = NULL, help = "second replicate path", metavar = "character"),
    make_option(c("-w", "--whitelist"), type="character", default = NULL, help = "gene filter, anything here considered as passing the filtering criteria", metavar = "character"),
    make_option(c("-l", "--aliasname"), type="character", default = NULL, help = "gff like file to convert coordinate to gene name", metavar = "character"),
    make_option(c("-o", "--output"), type="character", default = NULL, help = "output file", metavar = "character")
)
 
opt_parser = OptionParser(option_list= option_list)
opt = parse_args(opt_parser)
#End of opt section

compile1 <- read.table(opt$input1, header= TRUE, sep= "\t")
compile2 <- read.table(opt$input2, header= TRUE, sep= "\t")

#I did this because paternally imprinted genes are not required to be filtered
pcompile1 <- compile1[compile1$IMP == 'P',]
pcompile2 <- compile2[compile2$IMP == 'P',]

compileref1 <- read.table(opt$aliasname, header= FALSE, sep= "\t")
compileref2 <- read.table(opt$whitelist, header= TRUE, sep= "\t")
raw1 <- read.table(opt$raw1, header= TRUE, sep= "\t")
raw1 <- raw1[c('NAMES','EXON1','GENECAT1','PARSUM1','PARSUM2','TOTSUM1','TOTSUM2','PARFRAC1','PARFRAC2','IMP1','IMP2','genchitest1','genchitest2')]
raw2 <- read.table(opt$raw2, header= TRUE, sep= "\t")
raw2 <- raw2[c('NAMES','EXON1','GENECAT1','PARSUM1','PARSUM2','TOTSUM1','TOTSUM2','PARFRAC1','PARFRAC2','IMP1','IMP2','genchitest1','genchitest2')]
searchtable <- compileref1[match(compileref2$Names, compileref1[,1]),]

#Filtering the endosperm specific expressed genes
mcompile1 <- compile1[compile1$IMP == 'M',]
mcompile2 <- compile2[compile2$IMP == 'M',]
mcompile1 <- mcompile1[mcompile1$GENECAT %in% searchtable[,8],]
mcompile2 <- mcompile2[mcompile2$GENECAT %in% searchtable[,8],]

#write.table(mcompile1, file= "/home/rocky/Imprinting/mcompilefilr1.csv", sep= "\t", quote= FALSE, row.names= FALSE, col.names= TRUE)
#write.table(pcompile1, file= "/home/rocky/Imprinting/pcompilefilr1.csv", sep= "\t", quote= FALSE, row.names= FALSE, col.names= TRUE)
#write.table(mcompile2, file= "/home/rocky/Imprinting/mcompilefilr2.csv", sep= "\t", quote= FALSE, row.names= FALSE, col.names= TRUE)
#write.table(pcompile2, file= "/home/rocky/Imprinting/pcompilefilr2.csv", sep= "\t", quote= FALSE, row.names= FALSE, col.names= TRUE)

compile1 <- rbind(mcompile1,pcompile1)
compile2 <- rbind(mcompile2,pcompile2)

compiling <- rbind(compile1, compile2)
nrow(compiling)
compilingdup <- compiling[duplicated(compiling[,1:2]),]
compilingnodup <- compiling[!compiling$GENECAT %in% compilingdup$GENECAT,]
compiling <- compilingdup

#Please make notes that nodup suffixes here doesn't mean not duplicated. It means that the threshold has beed adjusted to accept that being imprinted in one library from each reciprocal crosses are enough justify imprinting.
compilingdup1 <- raw1[which(raw1$GENECAT1 %in% compilingdup$GENECAT),]
compilingdup2 <- raw2[which(raw2$GENECAT1 %in% compilingdup$GENECAT),]
compilingnodup1 <- raw1[which(raw1$GENECAT1 %in% compilingnodup$GENECAT),]
compilingnodup2 <- raw2[which(raw2$GENECAT1 %in% compilingnodup$GENECAT),]
compilingdup <- merge(compilingdup1, compilingdup2, by.x= "NAMES", by.y= "NAMES", all= TRUE)
compilingnodup <- merge(compilingnodup1, compilingnodup2, by.x= "NAMES", by.y= "NAMES", all= TRUE)
#write.table(compilingdup, file= "/home/rocky/Imprinting/compiledup.csv", sep= "\t", quote= FALSE, row.names= FALSE, col.names= TRUE)
#write.table(compilingnodup, file= "/home/rocky/Imprinting/compilenodup.csv", sep= "\t", quote= FALSE, row.names= FALSE, col.names= TRUE)
#write.table(compilingnodup1, file= "/home/rocky/Imprinting/compilenodup1.csv", sep= "\t", quote= FALSE, row.names= FALSE, col.names= TRUE)
#write.table(compilingnodup2, file= "/home/rocky/Imprinting/compilenodup2.csv", sep= "\t", quote= FALSE, row.names= FALSE, col.names= TRUE)
compilingnodupfil <- subset(compilingnodup, TOTSUM1.x >= 10 & TOTSUM2.x >= 10 & TOTSUM1.y >= 10 & TOTSUM2.y >= 10)
compilingnodupfil <- subset(compilingnodupfil, PARFRAC1.x >= 0.666667 & PARFRAC2.x >= 0.666667 & PARFRAC1.y >= 0.666667 & PARFRAC2.y >= 0.666667)
compilingnodupfil <- subset(compilingnodupfil, IMP1.x == IMP2.x & IMP1.y == IMP2.y & IMP1.x == IMP2.y)
#write.table(compilingnodupfil, file= "/media/diskb/rocky/rugraproj/paired_end/unique/compilenodupfilcrr2.csv", sep= "\t", quote= FALSE, row.names= FALSE, col.names= TRUE)

#Write duplicated genes in both replicate. nodup suffixes here means that these imprinted genes doesn't covered well in other replicates or termed as no evidence. No evidence as imprinted is not the same as proven to be not imprinted. It's still imprinted well in one replicate.
#Based on this, reported imprinted genes are combination of three files, dup, nodup1, and nodup2.
pcompilingdup <- compiling[compiling$IMP == "P",]
mcompilingdup <- compiling[compiling$IMP == "M",]

#This block is to add GeneID
compileref1 <- compileref1[,c(8,1)]
names(compileref1) <- c("GENECAT", "GeneID")
pcompilingdup <- merge(pcompilingdup, compileref1, by= "GENECAT", all.x= TRUE)
mcompilingdup <- merge(mcompilingdup, compileref1, by= "GENECAT", all.x= TRUE)

write.table(pcompilingdup, file= paste(opt$output,"p.csv", sep= ""), sep= "\t", quote= FALSE, row.names= FALSE, col.names= TRUE)
write.table(mcompilingdup, file= paste(opt$output,"m.csv", sep= ""), sep= "\t", quote= FALSE, row.names= FALSE, col.names= TRUE)

compilingnodup <- compiling[!compiling$GENECAT %in% compilingdup$GENECAT,]
compilingnodup1 <- compile1[which(compile1$GENECAT %in% compilingnodupfil$GENECAT1.x),]
compilingnodup2 <- compile2[which(compile2$GENECAT %in% compilingnodupfil$GENECAT1.y),]

pcompilingnodup1 <- compilingnodup1[compilingnodup1$IMP == "P",]
mcompilingnodup1 <- compilingnodup1[compilingnodup1$IMP == "M",]
pcompilingnodup2 <- compilingnodup2[compilingnodup2$IMP == "P",]
mcompilingnodup2 <- compilingnodup2[compilingnodup2$IMP == "M",]

#This block is to add GeneID
pcompilingnodup1 <- merge(pcompilingnodup1, compileref1, by= "GENECAT", all.x= TRUE)
mcompilingnodup1 <- merge(mcompilingnodup1, compileref1, by= "GENECAT", all.x= TRUE)
pcompilingnodup2 <- merge(pcompilingnodup2, compileref1, by= "GENECAT", all.x= TRUE)
mcompilingnodup2 <- merge(mcompilingnodup2, compileref1, by= "GENECAT", all.x= TRUE)

#write.table(pcompilingnodup1, file= "/media/diskb/rocky/rugraproj/paired_end/unique/pcompilenodupcrr2.csv", sep= "\t", quote= FALSE, row.names= FALSE, col.names= TRUE)
#write.table(mcompilingnodup1, file= "/media/diskb/rocky/rugraproj/paired_end/unique/mcompilenodupcrr2.csv", sep= "\t", quote= FALSE, row.names= FALSE, col.names= TRUE)
#write.table(pcompilingnodup2, file= "/media/diskb/rocky/rugraproj/paired_end/unique/pcompilenodupcrr2.csv", sep= "\t", quote= FALSE, row.names= FALSE, col.names= TRUE)
#write.table(mcompilingnodup2, file= "/media/diskb/rocky/rugraproj/paired_end/unique/mcompilenodupcrr2.csv", sep= "\t", quote= FALSE, row.names= FALSE, col.names= TRUE)
