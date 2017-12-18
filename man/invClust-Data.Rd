\name{geno}
\alias{geno}
\alias{annot}
\alias{roi}
\title{Demo data for invClust}

\description{
Data sets to illustrate the use of \code{invCLust}, \code{breed} and \code{invClustBreed}.
They correspond to the MKK subsample of the HapMapIII populations.  
}

\usage{
data(geno)
}

\format{
\code{geno} is a SnpMatrix of SNP genotypes, \code{annot} its SNP annotation data and \code{roi} the 
required imput for the genomic segment to be analized.
\code{haplos} is a data.frame of phased alleles, where columns are loci and rows encode subject chromsomes.
\code{bd} is a list and illustrated the output of \code{invClustBreed}
}

\examples{
data(geno)
geno
head(annot)
roi

}
