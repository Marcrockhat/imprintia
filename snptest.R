#MRH 18112013
#This script is used to identify maternal and paternal SNP and calculate the p-value on proportion test based on several methods.
rm(list = ls())
library('WRS')
snpdat <- read.table("/home/rocky/Imprinting/Cr48xCr75.csv", header= FALSE, sep= " ", fill= TRUE)
names(snpdat) = c("NAMES","POS", "NUMA", "NUMC", "NUMG", "NUMT", "NUMN")

#Revert the column names for technical reason, I'm using ATGC as my conventional way to name column
#flipping this order will change the outcome of the analysis since I'm using fixed index which related to the column name for the bases
#A is 1, T is 3, G is 9, and C is 27. If you want to change the order, please rewrite the index which is used to distinguish the maternal and paternal alleles
snpdat <- snpdat[c("NAMES","POS", "NUMA", "NUMT", "NUMG", "NUMC", "NUMN")]
snpdat[is.na(snpdat)] <- 0

#Extract additional information for each locus including parental allele type and gene name.
snpref <- read.table("/home/rocky/Imprinting/compfil48snp.csv", header= TRUE, sep= "\t", fill= TRUE)
snpdat$MAT <- snpref$MAT
snpdat$PAT <- snpref$PAT
genecat <- snpref$genecat

rm(snpref)
matsnp <- as.matrix(snpdat$MAT)
newmatsnp <- as.matrix(snpdat$MAT)
patsnp <- as.matrix(snpdat$PAT)
newpatsnp <- as.matrix(snpdat$PAT)
nallele <- 4

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
  matsnp <- cbind(matsnp, snpdat$MAT%%(nallele-1)^i)
}
for (i in 1:nallele){
  matsnp <- cbind(matsnp, matsnp[,i+1]%/%(nallele-1)^(i-1))
}

for (i in 1:nallele){
  patsnp <- cbind(patsnp, snpdat$PAT%%(nallele-1)^i)
}
for (i in 1:nallele){
  patsnp <- cbind(patsnp, patsnp[,i+1]%/%(nallele-1)^(i-1))
}

#This sets all the number to 1 we don't want to have 2 in the calculation since it will multiplied the number of alleles
matsnp <- cbind(matsnp,ceiling(matsnp[,6:9]/2))
patsnp <- cbind(patsnp,ceiling(patsnp[,6:9]/2))

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
snpdat$POS[snpdat$POS > 0] <- 1
snpdat <- cbind(snpdat, as.numeric(sumtype))
snpdat <- cbind(snpdat, as.numeric(newdmat$MATSUM))
snpdat <- cbind(snpdat, as.numeric(newdpat$PATSUM))
snpdat <- cbind(snpdat, as.numeric(sumpar))
names(snpdat) <- c("NAMES","POS", "NUMA", "NUMC", "NUMG", "NUMT", "NUMN", "GENECAT", "TYPE", "MATSUM", "PATSUM", "PARSUM")
snpdat <- snpdat[snpdat$POS > 0,]
snpdat <- snpdat[snpdat$TYPE > 0,]
snpdat <- snpdat[c("NAMES", "POS", "TYPE", "GENECAT", "MATSUM","PATSUM","PARSUM")]

matproptest <- NULL
patproptest <- NULL
parchitest <- NULL
matbintest <- NULL
mattwobintest <- NULL
pattwobintest <- NULL

summat <- snpdat$MATSUM
sumpat <- snpdat$PATSUM
sumpar <- snpdat$PARSUM

#Test the maternal binomial proportion test
for (i in 1:length(summat)){
  if(summat[i] > 0 & sumpar[i] > 0){
    mattest <- prop.test(summat[i], sumpar[i], p= (2/3))
    matproptest <- rbind(matproptest, mattest$p.value)
  } else {
    matproptest <- rbind(matproptest, NA)
  }
}

#Test the paternal binomial proportion test
for (i in 1:length(sumpat)){
  if(sumpat[i] > 0 & sumpar[i] > 0){
    pattest <- prop.test(sumpat[i], sumpar[i], p= (1/3))
    patproptest <- rbind(patproptest, pattest$p.value)
  } else {
    patproptest <- rbind(patproptest, NA)
  }
}

#Test the chisquare test
for (i in 1:length(sumpar)){
  if(summat[i] > 0 | sumpat[i] > 0){
    chitest <- chisq.test(c(summat[i],sumpat[i]), p= c(2/3, 1/3))
    parchitest <- rbind(parchitest, chitest$p.value)
  } else {
    parchitest <- rbind(parchitest, NA)
  }
}

#Test the maternal storer kim test
for (i in 1:length(sumpar)){
  if(summat[i] > 0 | sumpar[i] > 0){
    mattwotest <- twobinom(summat[i],sumpar[i], 2, 3)
    mattwobintest <- rbind(mattwobintest, mattwotest$p.value)
  } else {
    mattwobintest <- rbind(mattwobintest, NA)
  }
}

#Test the paternal storer kim test
for (i in 1:length(sumpar)){
  if(sumpat[i] > 0 | sumpar[i] > 0){
    pattwotest <- twobinom(sumpat[i],sumpar[i], 1, 3)
    pattwobintest <- rbind(pattwobintest, pattwotest$p.value)
  } else {
    pattwobintest <- rbind(pattwobintest, NA)
  }
}

snpdat <- cbind(snpdat, matproptest)
snpdat <- cbind(snpdat, patproptest)
snpdat <- cbind(snpdat, parchitest)
snpdat <- cbind(snpdat, mattwobintest)
snpdat <- cbind(snpdat, pattwobintest)

write.table(snpdat, file= "/home/rocky/Imprinting/Cr48xCr75calc.csv", sep= "\t", quote= FALSE, row.names= FALSE, col.names= TRUE)
