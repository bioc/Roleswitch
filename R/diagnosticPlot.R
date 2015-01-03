# diagnostic plot
# Input: a list object from roleswitch function

diagnosticPlot <- function(pred) {
  
  stopifnot(class(pred)=="ProMISe")
  
	par(mfrow=c(2,4))
	
	x.o <- pred$x.o
	z.o <- pred$z.o
	c <- pred$c
	x.t <- pred$x.t
	z.t <- pred$z.t
	p.x <- pred$p.x
	p.z <- pred$p.z
	p.xz <- pred$p.xz
	delta.p.all <- pred$delta.p.all
	
  i <- 1
  
	# obs mRNA expression
	color2D.matplot(x.o, extremes=c("white", "black"),
		main=sprintf("(%s) Obs. mRNA expr", LETTERS[i]), 
		axes=FALSE, xlab="", ylab="", show.values=2)
	
	axis(2,at=0.5:(nrow(x.o)-0.5),las=2,labels=nrow(x.o):1)
	
  i <- i + 1
  
	# obs miRNA expression
	color2D.matplot(z.o, extremes=c("white", "black"),
		main=sprintf("(%s) Obs. miRNA expr", LETTERS[i]), axes=FALSE, xlab="", ylab="",
		show.values=2)
		
	axis(2,at=0.5:(nrow(z.o)-0.5),las=2,labels=nrow(z.o):1)	
	
  i <- i + 1
  
	# obs seed matrix
	color2D.matplot(c, extremes=c("white", "black"),
		main=sprintf("(%s) Obs. seed match", LETTERS[i]), axes=FALSE, xlab="miRNA", ylab="mRNA", show.values=2)
	
	axis(1,at=0.5:(nrow(z.o)-0.5),las=1,labels=1:nrow(z.o))
	
	axis(2,at=0.5:(nrow(x.o)-0.5),las=2,labels=nrow(x.o):1)
	
  i <- i + 1
  
	# inferred mRNA expression
	color2D.matplot(x.t, extremes=c("white", "black"),
		main=sprintf("(%s) Est. mRNA total", LETTERS[i]), axes=FALSE, xlab="", ylab="",
		show.values=2)
	
	axis(2,at=0.5:(nrow(x.o)-0.5),las=2,labels=nrow(x.o):1)
		
  i <- i + 1
  
	# inferred miRNA-mRNA Probability Matrix
	color2D.matplot(p.x, extremes=c("white", "black"),
		main=sprintf("(%s) mRNA competition", LETTERS[i]), axes=FALSE, xlab="miRNA", ylab="mRNA",
		show.values=2)
	
	axis(1,at=0.5:(nrow(z.o)-0.5),las=1,labels=1:nrow(z.o))
	
	axis(2,at=0.5:(nrow(x.o)-0.5),las=2,labels=nrow(x.o):1)
	
  i <- i + 1
  
	# inferred mRNA-miRNA Probability Matrix
	color2D.matplot(p.z, extremes=c("white", "black"),
		main=sprintf("(%s) miRNA competition", LETTERS[i]), axes=FALSE, xlab="miRNA", ylab="mRNA",
		show.values=2)
	
	axis(1,at=0.5:(nrow(z.o)-0.5),las=1,labels=1:nrow(z.o))
	
	axis(2,at=0.5:(nrow(x.o)-0.5),las=2,labels=nrow(x.o):1)
	
  i <- i + 1
  
	# inferred final Probability Matrix
	color2D.matplot(p.xz, extremes=c("white", "black"),
		main=sprintf("(%s) Joint competition", LETTERS[i]), axes=FALSE, xlab="miRNA", ylab="mRNA",
		show.values=2)
	
	axis(1,at=0.5:(nrow(z.o)-0.5),las=1,labels=1:nrow(z.o))
	
	axis(2,at=0.5:(nrow(x.o)-0.5),las=2,labels=nrow(x.o):1)

  i <- i + 1
  
	# convergence progress
	plot(delta.p.all[-1], ylab="max(abs(p.x - p.x.prev))", 
       xlab="iteration", main=sprintf("(%s) Convergence", LETTERS[i]))
}



