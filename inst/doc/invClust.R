### R code from vignette source 'invClust.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: loadpck
###################################################
library(invClust)


###################################################
### code chunk number 2: SNPstats (eval = FALSE)
###################################################
## source("http://bioconductor.org/biocLite.R")
## biocLite("snpStats")


###################################################
### code chunk number 3: load
###################################################
library(snpStats)
path<-system.file("extdata", package = "invClust")
genofile<-file.path(path,"genosPLINK")
#genofile is the path where the genoPLINK demo files are stored
#use your own path for your own data 
geno.data<-read.plink(genofile)


###################################################
### code chunk number 4: geno
###################################################
geno<-geno.data$genotypes
geno


###################################################
### code chunk number 5: anot
###################################################
annot.read<-geno.data$map
annot<-annot.read[,c(1,2,4)]
head(annot)


###################################################
### code chunk number 6: check
###################################################
identical(annot[,2],colnames(geno))


###################################################
### code chunk number 7: check
###################################################
geno.mat<-matrix(c(0,0,1,2,0,0),ncol=2)
rownames(geno.mat)<-c("sub1","sub2","sub3")
colnames(geno.mat)<-c("rs1","rs2")
geno.raw<-matrix(as.raw(geno.mat+1),ncol=ncol(geno.mat))
geno.new<-new("SnpMatrix", geno.mat)
geno.new


###################################################
### code chunk number 8: check
###################################################
roi<-data.frame(chr=8,LBP=7934925, RBP=11824441, reg= "inv1")


###################################################
### code chunk number 9: call (eval = FALSE)
###################################################
## invcall<-invClust(roi=roi, wh = 1, geno=geno, annot=annot, dim=2)


###################################################
### code chunk number 10: check (eval = FALSE)
###################################################
## invcall<-invClust(roi="ROI.txt", wh = 1, geno=geno, annot=annot, dim=2)


###################################################
### code chunk number 11: calltrue
###################################################
quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}
invcall<-quiet(invClust(roi=roi, wh = 1, geno=geno, annot=annot, dim=2))


###################################################
### code chunk number 12: print
###################################################
invcall


###################################################
### code chunk number 13: plot
###################################################
invcall
plot(invcall)


###################################################
### code chunk number 14: invgeno
###################################################
inv<-invGenotypes(invcall)
head(inv)


###################################################
### code chunk number 15: invUnc
###################################################
invUnc<-invcall["genotypes"]
head(invUnc)


###################################################
### code chunk number 16: snpassoc (eval = FALSE)
###################################################
## install.packages("SNPassoc")


###################################################
### code chunk number 17: snpassoc
###################################################
library(SNPassoc)


###################################################
### code chunk number 18: BMI
###################################################
path<-system.file("extdata", package = "invClust")
phenofile<-file.path(path,"BMI.txt")
BMI<-read.delim(phenofile,as.is=TRUE)
head(BMI)


###################################################
### code chunk number 19: BMI
###################################################
identical(BMI$ID, names(inv))
data<-cbind(BMI,inv)


###################################################
### code chunk number 20: setup
###################################################
data.end<-setupSNP(data,colSNPs=3)
head(data.end)
class(data.end$inv)


###################################################
### code chunk number 21: setup
###################################################
association(BMI~inv,data.end)


