\name{tcga_gbm_testdata}
\alias{tcga_gbm_testdata}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
Test data of miRNA and mRNA expression vector.
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
Test data of miRNA and mRNA expression vector from the same human individial that is used in the vignette.
}
\usage{data(tcga_gbm_testdata)}
\format{
  A data frame with 11884 mRNA expression and 373 miRNA expression  on the following 2 variables.
  \describe{
    \item{\code{x}}{a numeric vector for expression of 11884 mRNAs}
    \item{\code{y}}{a numeric vector for expression of 373 miRNA}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
The miRNA and mRNA expression data for the same individual (barcode ID: TCGA-02-0001-01) were downloaded from TCGA (The Cancer Genome Altas) GBM (Glioblastoma multiforme). To eliminate negative values in the expression matrices, we linearly transformed all of the data to positive scale with range between 0 and the maximum of the positive values using \code{\link[scales]{rescale}} R function. The resulting data were further processed by filtering out miRNAs or mRNAs without any seed/seed match based on the Microcosm database. As a result, the test data contain the expression for 11884 distinct mRNA and distinct 373 miRNA.
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
\url{https://tcga-data.nci.nih.gov}
}
\references{
%%  ~~ possibly secondary sources and usages ~~
Griffiths-Jones S, Saini HK, van Dongen S, Enright AJ (2008) miRBase: tools for microRNA genomics. Nucleic acids research 36: D154-8
}
\examples{
data(tcga_gbm_testdata)
}
\keyword{datasets}
