\name{diagnosticPlot}
\alias{diagnosticPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create diagnostic plot for understanding the Roleswitch model outputs.
}
\description{
Create diagnostic plot for understanding the \code{\link{roleswitch}} outputs. Create a 2 by 4 panels of plots. From left to right, the top panel displays the observed N mRNA and M miRNA expression (\code{x.o} and \code{z.o}), the N by M seed-match matrix (\code{c}), and the inferred total mRNA expression; the bottom panel displays the inferred probability of the M miRNAs targeting the N mRNA (miRNA-mRNA; \code{p.x}), the probability of the N mRNA "targeting" the M miRNAs (mRNA-miRNA; \code{p.z}), the dot product of the above two matrices (Joint) and the convergence rate (\code{delta.p.all}). 
}
\usage{
diagnosticPlot(pred)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pred}{
Results obtained from \code{\link{roleswitch}}.
}
}
\references{
Li, Y., ..., Zhang, Z., Inference of personalized miRNA-mRNA interactions toward redefining cancer signatures (in preparation).
}
\author{
Yue Li
}

\seealso{
\code{\link{roleswitch}}
}
\examples{
x.o <- matrix(abs(rnorm(10, mean=3)), 
  dimnames=list(c(1:10),"mRNA")) # mRNA expression
  
z.o <- matrix(abs(rnorm(4, mean=3)),
  dimnames=list(c(1:4),"miRNA"))  # miRNA expression
  
c <- matrix(rpois(40, lambda=3),nrow=nrow(x.o), 
  dimnames=list(c(1:10),c(1:4))) # seed match matrix

rs.pred <- roleswitch(x.o, z.o, c)

diagnosticPlot(rs.pred)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{plot}
