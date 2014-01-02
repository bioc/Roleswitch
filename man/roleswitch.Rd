\name{roleswitch}
\alias{roleswitch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Infer miRNA-mRNA interactions using paired expression data from a single sample.
}
\description{
Infer Probabilities of MiRNA-mRNA Interaction Signature (ProMISe) using paired expression data from a single sample. Roleswitch operates in two phases by inferring the probabilities of mRNA (miRNA) being the targets ("targets") of miRNA (mRNA), taking into account the expression of all of the mRNAs (miRNAs) due to their potential competition for the same miRNA (mRNA). Due to mRNA transcription and miRNA repression events simultaneously happening in the cell, Roleswitch assumes that the total transcribed mRNA levels are higher than the observed (equilibrium) mRNA levels and iteratively updates the total transcription of each mRNA targets based on the above inference. NB: in the paper, we used ProMISe as both the model name and inferred score name.
}
\usage{
roleswitch(x.o, z.o, c, maxiter = 200, tol = 1e-05,
  eta.z = 0.001, expected.total = 1.3, verbose = TRUE, 
  annotation.db, probe2genesymbol=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x.o}{
A numeric vector as the observed expression of N mRNAs. NOTE: rownames is required for x.o to map each mRNA to the rownames of the seed match matrix.
}
  \item{z.o}{
A numeric vector as the observed expression of M miRNAs. NOTE: rownames is required for z.o to map each miRNA to the colnames of the seed match matrix.
}
  \item{c}{
A numeric N x M matrix of integers representing the seed match matrix between N mRNA and M miRNA. NOTE: dimnames is required for the seed match matrix to map each mRNA and each miRNA to its rownames and colnames, respectively.
}
  \item{maxiter}{
The maximum number of iterations before terminating roleswitch inference (default: 200).
}
  \item{tol}{
The threshold on the largest absolute difference between the current and the previous probabilities of miRNA-mRNA interactions (default: 1e-5).
}
  \item{eta.z}{
A scalar decimal value specifying the amount of update applied to the total transribed mRNA, which is inferred durin the iteration (default: 1e-3).
}
  \item{expected.total}{
The ratio of total transcription over observed transcription (default: 1.3).
}
  \item{verbose}{
Display progress at each iteration.
}
  \item{annotation.db}{
Charater string specifying the name of the annotation package for microarray platform if the input \code{x.o} is an \code{\link[Biobase]{eSet}} or \code{\link[Biobase]{ExpressionSet}} object. This is optional only if annotation slot in the eSet \code{x.o} is defined.
}
  \item{probe2genesymbol}{
Whether to convert probe id to gene symbol (Default: TRUE). This only applies when \code{x.o} is an eSet/ExpressionSet. Probe ID are usually the \code{featureNames}, an attribute under the class \code{\link[Biobase]{eSet}}. If TRUE, conversion from probe id to gene symbol is performed automatically.
}
  \item{\dots}{
Arguments passed to \code{\link{getSeedMatrix}}.
}
}
\details{
The model assumes total expression of mRNA is unobserved and higher than the observed corresponidng expression due to RNA degradation induced by miRNA-mRNA interaction. The general algorithm is outlined as follows:

(1) Infer mRNA i targeted by miRNA k taking into account the hidden total expression of 1...N mRNA and miRNA k

(2) Estimate total transcription level of mRNA i

(3) Infer miRNA k "targeted" by mRNA i taking into account 1...M miRNA and mRNA i expression

(4) Repeat 1-3 until convergence

User provide roleswitch an N x M seed-match matrix containg the number of target sites for each mRNA i and miRNA k. Otherwise, \code{\link{getSeedMatrix}} will be used to retrieve seed-match matrix. The output ProMISe is one of the matrices: mRNA competition; miRNA competition; joint competition. We recommand using mRNA competition or joint competition.
}
\value{
An object list (defined class: ProMISe) containing the following items:
\item{x.t:}{infered total N x 1 mRNA expression vector}
\item{p.x:}{infered N x M miRNA-mRNA probability matrix (mRNA competition)}
\item{p.z:}{infered N x M mRNA-miRNA probability matrix (miRNA competition)}
\item{p.xz:}{Dot prodct of p.x and p.z (joint competition).}
\item{x.o, z.t, c:}{The same as the inputs.}
\item{delta.p.all}{Difference between p.x at two adjacent iterations to monitor the progress of convergence.}
}
\references{
Li, Y., ..., Zhang, Z., Infer probabilistic miRNA-mRNA interaction signatures in cancers: a role-switch approach (in preparation).
}
\author{
Yue Li
}
\note{
Warning may be issued if the input seed matrix contain miRNA or mRNA that have zero seed or seed match for any mRNA or miRNA, respectively. Nonetheless, the outputs ProMISe will conform the original input mRNA and miRNA vector in matching their corresponding IDs.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{getSeedMatrix}}
}
\examples{
x.o <- matrix(abs(rnorm(10, mean=3)), 
  dimnames=list(c(1:10),"mRNA")) # mRNA expression
  
z.o <- matrix(abs(rnorm(4, mean=3)),
  dimnames=list(c(1:4),"miRNA"))  # miRNA expression
  
c <- matrix(rpois(40, lambda=3),nrow=nrow(x.o), 
  dimnames=list(c(1:10),c(1:4))) # seed match matrix

rs.pred <- roleswitch(x.o, z.o, c)

rs.pred$p.xz
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{roleswitch}
\keyword{miRNA}
\keyword{microRNA}
\keyword{expression}


