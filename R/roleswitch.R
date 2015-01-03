# Roleswitch model in miRNA target prediction
# The model assumes total expression of mRNA is unobserved and higher than the observed corresponidng expression due to RNA degradation induced by miRNA-mRNA interaction
# (1) Infer mRNA i targeted by miRNA k taking into account the hidden total expression of 1...N mRNA and miRNA k
# (2) Estimate total transcription level of mRNA i
# (3) Infer miRNA k "targeted" by mRNA i taking into account 1...M miRNA and mRNA i expression
# (4) Repeat 1-3 until convergence

# Input:
# x.o: observed N x 1 mRNA expression vector
# z.o: observed M x 1 miRNA expression vector
# c: observed N x M seed match matrix
# cellcap: total expression w/o miRNA regulatory effects (default: 1.3 * x.o)
# tol: theshold for probability matrix convergence
# roleswitch: whether to estimate the total miRNA abundance or assume the observed miRNA abundance is the total miRNA abudance
# ...: parameters passed to getSeedMatrix

# Output:
# ProMISe: A list containing
# x.t: infered total N x 1 mRNA expression vector
# z.t: infered total M x 1 miRNA expression vector
# p.x: infered N x M miRNA-mRNA probability matrix (mRNA plays as targets of miRNA)
# p.z: infered N x M mRNA-miRNA probability matrix (miRNA plays as targets of mRNA)
# p.xz: element-wise product of p.x and p.z (final inference)

roleswitch <- function(x.o, z.o, c, maxiter=200, tol=1e-5, 
	eta.z = 0.001, expected.total=1.3, # normalize=FALSE,
	verbose=TRUE, annotation.db, probe2genesymbol=TRUE, ...) {	
  
  # process eSet/ExpressionSet
  if(class(x.o)=="eSet" || class(x.o)=="ExpressionSet") {
    
    genes <- featureNames(x.o)        
    
    if(probe2genesymbol) {
      
      if(missing(annotation.db) & is.null(annotation(x.o)))
        stop("annotation.db or annotation is needed to map probe to gene id")
      
      if(missing(annotation.db) & !is.null(annotation(x.o)))
        annotation.db <- sprintf("%s.db", annotation(x.o))
                   
      if(is.null(featureNames(x.o)))
        stop("featureNames is missing in the eSet/ExpressionSet to map probe id to gene id")            
      
      require(annotation.db, character.only = TRUE) ||
        stop(sprintf("%s package must be installed", annotation.db))
      
      probeID <- genes      
      
      con_comm <- paste(as.character(strsplit(annotation.db, split = ".db")), 
                        "dbconn", sep = "_")
      
      usedDB <- eval(as.name(con_comm))
      
      con <- usedDB()
      
      annot.s <- dbReadTable(con, "gene_info") #[, c("X_id", "symbol")]
      
      annot.e <- dbReadTable(con, "genes")
      
      annot <- merge(merge(annot.s, annot.e), dbReadTable(con, "probes"))
      
      genes <- annot[match(probeID, annot$probe_id),"symbol"]      
    }      
    
    if(ncol(x.o) > 1) warning("mRNA expression matrix has >1 column. Only first column is used.")
    
    x.o <- exprs(x.o)[,1,drop=F]       
        
    if(any(duplicated(genes))) {
      
      message("Average probe values were taken over the same gene.")
      
      df <- aggregate(data.frame(expr=as.numeric(x.o)), list(genes=genes), mean)
      
      x.o <- as.matrix(df$expr)
      
      rownames(x.o) <- as.matrix(df$genes)
      
    } else {      
      rownames(x.o) <- as.matrix(genes)      
    }
    
    x.o <- x.o[!is.na(rownames(x.o)),,drop=F]
  }
  
  if(ncol(z.o) > 1) {
    warning("miRNA expression matrix has >1 column. Only first column is used.")
    z.o <- z.o[,1,drop=F]
  }
  
  if(missing(c)) {
    
    c <- getSeedMatrix(convert2genesymbol=TRUE, ...)
        
  } else {
    
    stopifnot(all(dim(c)==c(nrow(x.o),nrow(z.o))))
    
    stopifnot(all(sort(rownames(c))==sort(rownames(x.o))))
    
    stopifnot(all(sort(colnames(c))==sort(rownames(z.o))))    
  }
  
  stopifnot(!is.null(rownames(x.o)) & !is.null(rownames(z.o)) & !is.null(dimnames(c)))  
  
  if(ncol(x.o) > 1) warning("mRNA expression matrix has >1 column. Only first column is used.")
  
  x.o <- x.o[,1,drop=F]
  
  if(any(x.o<0)) {
    
    message("Negative values detected in mRNA expression are set to 0")
    
    x.o[x.o<0] <- 0    
  }
  
  if(any(z.o<0)) {    
    
    message("Negative values detected in miRNA expression are set to 0")
    
    z.o[z.o<0] <- 0    
  }
  
  if(any(duplicated(rownames(z.o)))) {
    
    message("Average values were taken over the same miRNA.")
    
    df <- aggregate(data.frame(expr=as.numeric(z.o)), list(mirna=rownames(z.o)), mean)
    
    z.o <- as.matrix(df$expr)
    
    rownames(z.o) <- as.matrix(df$mirna)
  }
  
  x.o0 <- x.o
  z.o0 <- z.o  
  
	genes <- rownames(x.o)
	mirna <- rownames(z.o)
	
	c0 <- c[match(rownames(x.o), rownames(c)), 
		match(rownames(z.o), colnames(c))]
  
  dimnames(c0) <- list(rownames(x.o), rownames(z.o))    
			
	c0[is.na(c0)] <- 0
  
	c <- c0
	
	# remove mRNA (miRNA) that have no targeting miRNA (mRNA targets)
	rowidx <- which(apply(c, 1, sum) > 0)
	
	colidx <- which(apply(c, 2, sum) > 0)
	
	c <- c[rowidx,colidx]
	
	x.o <- x.o[rowidx,,drop=F]
	z.o <- z.o[colidx,,drop=F]
		
	### Begin role-switch estimation ###
	# initialize x.t and z.t by x.o and z.o
	x.t <- x.o
	z.t <- z.o
		
	x.total <- expected.total * sum(x.o) # total amount of mRNA (constant)
		
	t <- 1	
	converged <- FALSE
	
	p.x.prev <- matrix(0, nrow(x.o), nrow(z.o))	
	p.z.prev <- matrix(0, nrow(x.o), nrow(z.o))
	
	delta.p.all <- 1
	
	if(verbose) message(sprintf("\nStart roleswitch with %s miRNA and %s mRNA", 
		nrow(z.o), nrow(x.o)))
	
	while(t < maxiter & !converged) {
										
		t <- t + 1
		
		### miRNA competition for mRNA ###
		# total _e_xpressed _s_ites per miRNA
		ex.total <- matrix(rep(t(c) %*% x.t, nrow(x.t)), nrow(x.t), byrow=T)
		
		# expressed sites
		ex <- c * matrix(rep(x.t, nrow(z.t)), ncol=nrow(z.t))
			
		# remaining expressed sites = total expressed sites - expressed sites
		rex <- ex.total - ex
				
		# prob of mRNA target of miRNA k given total mRNA and miRNA k total
		# p(t_i|X^t,z_k^t)
		p.x <- 1 - exp(matrix(rep(z.t, nrow(x.t)), ncol=nrow(z.t), byrow=T) * 
					(log(rex) - log(ex.total)))
					
		p.x[is.infinite(p.x) | is.nan(p.x)] <- 0
			
		# delta(x_i|k): amount of mRNA i reduced by miRNA k
		x.d <- matrix(rep(x.t,nrow(z.t)),ncol=nrow(z.t)) * p.x * eta.z
						
		# update mRNA total
		x.t <- x.o + apply(x.d, 1, sum)
		
		# normalize by total transcribed (limited by cellular capacity)
		x.t <- (x.t * x.total)/sum(x.t)
		
		### mRNA competition for miRNA ###		
		# total expressed miRNA per mRNA (ei: _e_xpression m_i_RNA)
		ez.total <- matrix(rep(c %*% z.t, nrow(z.t)), ncol=nrow(z.t))
		
		# expressed miRNA k per mRNA
		ez <- c * matrix(rep(t(z.t),nrow(x.t)), nrow=nrow(x.t), byrow=T)
		
		# remaining expressed miRNA for mRNA to "target" =
		# total expressed miRNA per mRNA - expressed miRNA k per mRNA
		rez <- ez.total - ez	
			
		# prob of miRNA k given mRNA i and total miRNA
		# p(t_i|X^t,z_k^t)
		p.z <- 1 - exp(matrix(rep(x.t, nrow(z.t)), ncol=nrow(z.t)) * 
					(log(rez) - log(ez.total)))
					
		p.z[is.infinite(p.z) | is.nan(p.z)] <- 0									
		
		if(verbose & t > 1) message(sprintf("%d: max(p.x-p.x.prev)=%.5f", t-1, max(abs(p.x - p.x.prev))))
		
		delta.p <- max(abs(p.x - p.x.prev))
		
		delta.p2 <- max(abs(p.z - p.z.prev))
		
		delta.p.all <- c(delta.p.all, delta.p)
				
		if(delta.p < tol & delta.p2 < tol) {			
			
			converged <- T
			
		} else {
														
			p.x.prev <- p.x
			p.z.prev <- p.z
		}
	}
	
	# recover original list
	if(!identical(rownames(p.x), genes) || !identical(colnames(p.x), mirna)) {
		
		if(verbose) message("Some genes or miRNA are left out in calculation\nb/c they have zero target sites or targets!\nTheir probabilities are set to zero in the output matrices")
		
		p.x <- p.x[match(genes, rownames(p.x)), match(mirna, colnames(p.x))]
		p.z <- p.z[match(genes, rownames(p.z)), match(mirna, colnames(p.z))]
				
		p.x[is.na(p.x)] <- 0
		p.z[is.na(p.z)] <- 0
				
		dimnames(p.x) <- dimnames(c0)
		dimnames(p.z) <- dimnames(c0)
		
		x.t <- x.t[match(genes, rownames(x.t)),,drop=F]
		
		x.t[is.na(x.t)] <- x.o[is.na(x.t)]
		
		rownames(x.t) <- genes
	}
	
	p.xz <- p.x*p.z
	
	# if(normalize) {
		# p.xz <- normalize.p(p.xz, 2)
		# p.x <- normalize.p(p.x, 2)
		# p.z <- normalize.p(p.z, 2)		
	# }		
	
	promise <- list(x.t=x.t, z.t=z.o0, p.x=p.x, p.z=p.z, p.xz=p.xz, 
       c=c0, x.o=x.o0, z.o=z.o0, delta.p.all=delta.p.all[-1])
  
  class(promise) <- "ProMISe"
  
  promise
}

# normalize.p <- function(p, mydim) {
	
	# p[is.nan(p)] <- 0
	
	# s <- apply(p,mydim,sum)
	
	# if(mydim==1) {
				
		# p <- p/repmat(as.matrix(s),1,ncol(p))
		
		# # p[is.nan(p)] <- 1/ncol(p)
		
	# } else if(mydim==2) {
				
		# p <- p/repmat(s,nrow(p),1)
		
		# # p[is.nan(p)] <- 1/nrow(p)
	# }		
	
	# p[is.nan(p)] <- 0
	
	# p
# }


