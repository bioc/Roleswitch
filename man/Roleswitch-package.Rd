\name{Roleswitch-package}
\alias{Roleswitch-package}
\alias{Roleswitch}
\docType{package}
\title{
Infer miRNA-mRNA interactions using paired expression data from a single sample
}
\description{
Infer Probabilities of MiRNA-mRNA Interaction Signature (ProMISe) using paired expression data from a single sample. Roleswitch operates in two phases by inferring the probabilities of mRNA (miRNA) being the targets ("targets") of miRNA (mRNA), taking into account the expression of all of the mRNAs (miRNAs) due to their potential competition for the same miRNA (mRNA). Due to mRNA transcription and miRNA repression events simultaneously happening in the cell, Roleswitch assumes that the total transcribed mRNA levels are higher than the observed (equilibrium) mRNA levels and iteratively updates the total transcription of each mRNA targets based on the above inference. NB: in the paper, we used ProMISe as both the model name and inferred score name.
}
\details{
\tabular{ll}{
Package: \tab Roleswitch\cr
Type: \tab Package\cr
Version: \tab 1.1.2\cr
Date: \tab 2013-12-20\cr
License: \tab GPL-2\cr
}
The main function \code{\link{roleswitch}} takes as inputs the mRNA and miRNA expression as numerical vector and seedmatch (integer) matrix provided by the user or generated from \code{\link{getSeedMatrix}}. The function then outputs a list with p.xz containing the ProMISe as the final results and other results for diagnostic purpose.
}
\author{
Yue Li
Maintainer: Yue Li <yueli@cs.toronto.edu>
}
\references{
Li, Y., ..., Zhang, Z., Inference of personalized miRNA-mRNA interactions toward redefining cancer signatures (in preparation).

Arvey, A., Larsson, E., Sander, C., Leslie, C. S., & Marks, D. S. (2010). Target mRNA abundance dilutes microRNA and siRNA activity. Molecular systems biology, 6, 1-7. doi:10.1038/msb.2010.24
}

\keyword{Roleswitch}

\seealso{
\code{\link{roleswitch}}
}

\examples{
library(Roleswitch)
ls("package:Roleswitch")
}




