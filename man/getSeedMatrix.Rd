\name{getSeedMatrix}
\alias{getSeedMatrix}
\title{
Get seed-match matrix between defined mRNA and miRNA in an organism.
}
\description{
Given N mRNA and M miRNA IDs or simply the species common names (e.g., human), obtain the N by M seed match matrix as a the number of target sites each mRNA has for each miRNA.
}
\usage{
getSeedMatrix(mRNA, miRNA, species = "human", 
  id_type = "ensembl_transcript_id", mRNA_id_type = id_type,
  miRNA_id_type = id_type, longest3utr = TRUE, 
  biomart = "ensembl", dataset = "hsapiens_gene_ensembl",
  returnGeneInfo = FALSE, convert2genesymbol = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mRNA}{
A character vector of N mRNA ids.
}
  \item{miRNA}{
A character vector of M miRNA ids
}
  \item{species}{
Common names for a species. Currently only human and mouse are supported as precompiled target site information.
}
  \item{id_type}{
A string specifying the id type used for both the mRNAs and miRNA.
}
  \item{mRNA_id_type}{
A string specifying the id type used for the mRNAs.
}
  \item{miRNA_id_type}{
A string specifying the id type used for the miRNAs.
}
  \item{longest3utr}{
For genes having multiple transcripts, whetehr to use only the transcript with the longest 3'UTR (default: TRUE).
}
  \item{biomart}{
Database for biomart, which is ued to obtain transcript information using \code{\link[biomaRt]{getBM}} (default: ensembl).
}
  \item{dataset}{
Dataset used to query the biomart database using \code{\link[biomaRt]{getBM}} (default: hsapiens_gene_ensembl).
}
  \item{returnGeneInfo}{
Binary indicator to return gene information besides seed matrix (default: FALSE); if TRUE, then a list containing seed matrix and gene info is return; otherwise just the seed matrix.
}
  \item{convert2genesymbol}{
Whether to convert id such as ensembl_gene_id to gene symbols as row names of the seed matrix.
}
  \item{\dots}{
Paramters passed to \code{\link[biomaRt]{getBM}}.
}
}
\details{
Retrieve and process target site information to generate a N by M matrix representing the number of target sites of mRNA i for miRNA k. If species is specified, then the suggested data package RoleswitchData will be loaded and the pre-compiled seed matrix is used. Currently, only human and mouse are supported with this option. Otherwise, download the sequences based on the specified mRNA and miRNA IDs and obtain the seed matches using \code{\link[microRNA]{seedRegions}}.

}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{seed match matrix }{numeric matrix containing the number of target sites for each miRNA and mRNA pairs}
\item{gene info}{a data.frame containing miRNA id type, ensembl gene id, gene symbol, start and end of 3'UTR (only returned when \code{returnGeneInfo} is TRUE)}
}
\references{
miRBase: tools for microRNA genomics. (2008). miRBase: tools for microRNA genomics., 36(Database issue), D154-8. doi:10.1093/nar/gkm952

R. Gentleman and S. Falcon (2013). microRNA: Data and functions for dealing with microRNAs. R package version 1.18.0.
}
\author{
Yue Li
}
\note{
This is just a convenience function. Users are encouraged to construct the most up-to-date seed match matrix on their own from other source without using this function.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{roleswitch}},\code{\link[biomaRt]{getBM}},\code{\link[microRNA]{seedRegions}}
}
\examples{

seedMatrix.human <- getSeedMatrix()

head(seedMatrix.human)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ miRNA }
\keyword{ seedmatch }% __ONLY ONE__ keyword per line
