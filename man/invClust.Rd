\name{invClust}
\alias{invClust}
\alias{plot.invClust}
\alias{"[.invClust"}
\alias{getGenotypes}
\alias{invQuality }
\alias{print.invClust}

\title{
\code{invClust} Class and Methods
}
\description{
Mixture model fitting with Hardy Weinberg Equilibrium and population stratification to infer
haplotype (inversion) alleles.   
}
\usage{
invClust(roi, wh = 1, geno, annot, SNPtagg="n", SNPsel=1:ncol(geno), method=1, dim = 1, pc = 0, ngroups = 1, ...)
}
\arguments{     
  \item{roi}{text file of data.frame with Region of Interest information. Four columns 
  are required: chr, LBP, RBP, reg, with chromosome, left break point, right break point and a 
  character "reg" that identifies the inversion}
  \item{wh}{which ROI (row of roi) to be considered for computation.}     
  \item{geno}{genotypes in \code{snpMatrix} format of \code{snpStats}}
  \item{annot}{snp annotation as .map PLINK format, or \code{read.plink()$map} from snpStats. The columns
  "chromosome", "names", "position" are required and have to conform to this naming.}
  \item{SNPtagg}{set SNPtagg="y" to use tagg SNPs in roi for tagging the haplotype groups}
  \item{SNPsel}{vector with snps to be selected for computation}
  \item{method}{method=1 performs EM algorithm for three genotypes, method=2, performs a clustering within the 
   the genotypes for an additional third haplotype}   
  \item{dim}{ either 1 or 2 indicating the number fo MDS compoents to be used}
  \item{pc}{ if population clustering is to be performed, first component of a genome-wide PCA of geno}
  \item{ngroups}{maximum number of subpopualtions to be considered}
  \item{\dots}{control arguments for the EM algorithm: \code{it} (\code{it=1000}) maximum 
  number fo iterations, \code{tol} (\code{tol}=10e-5) convergence tolerance}
}
\details{
\code{invClust} computes the biallelic haplotypes in Hardy Weinberg Equilibrium (with the possibility
of clustering by geographical subpopulations) that may underlie an inversion event. It fits a mixture 
model with an expectation maximization routine, only controled by convergence criteria. Initial conditions
are general for a wide range of cases. Clustering is performed in 1 or 2 dimensions of mutidimeansional
scaling (argument \code{dim}) and, if geographical subpopulation is considered, the first 
component of a genome wide PCA (argument \code{pc}). In this last case, a visualization of the 
PCA analysis can inform on the suitable number of groups to be considered (argument \code{ngroup}). 

Each subject in the sample is assigned a probability to a given genotype (NN, NI, II) which 
can be recovered by \code{x["genotypes"]}, where \code{x} is of class \code{invClust} ( e.g.
the result of an \code{invClust} call). Most probable genotypes can be extracted with \code{getGenotypes(x)}. 
In a similar way, if subpopulation classification is considered,
probability for group membership is recovered with \code{x["groups"]}. 

Plots are also implemented for this class, 
\code{plot(x)} will display the clustered data on the fitted distribution, according to 
dimensions used. For inclussion of subpopulation classification, selection of the marginals 
can be done though a plot argument \code{wh=c("yy","xy")}. \code{wh="like"} plots the likelihood
with respect to the number of itereations only for \code{dim=1}.  

A useful quantity is a quality score (\code{getQuality(x)}) that computes the overlap integral of the 
cluster components, a value of 1 gives no overlap while 0 refers to complete overlap. 

}
\value{

\item{EMestimate }{ List with fitted parameters}
\item{datin }{ List with data used to fit the model: \code{x} for firt PCA component
and \code{y} MDS components }

}

\author{
Alejandro Caceres}
\examples{
data(geno)
inv<-invClust(roi=roi,wh=1,geno=geno,annot=annot,dim=1)
plot(inv)
head(inv["genotypes"])
}
\keyword{ invClust }
\keyword{ inversion }
