# Given mRNA and miRNA ID, count the number of seed matches
# for each mRNA-miRNA pair

getSeedMatrix <- function(mRNA, miRNA, species="human", id_type="ensembl_transcript_id", mRNA_id_type=id_type, miRNA_id_type=id_type, longest3utr=TRUE, biomart="ensembl", dataset="hsapiens_gene_ensembl", returnGeneInfo=FALSE, convert2genesymbol=TRUE, ...) {
	    
	if(species=="human") {  
    
    data(hsTargets, envir = environment())
    
    hsTargets <- get("hsTargets", envir = environment())
    
	seedMatrix <- with(hsTargets, table(target, name))
		
	} else {
	
		if(species=="mouse") {
      
		  data(mmTargets, envir = environment())
      
		  mmTargets <- get("mmTargets", envir  = environment())
		  
		  seedMatrix <- with(mmTargets, table(target, name))
			
		} else {
      
			mart <- useMart(biomart=biomart, dataset=dataset, ...)
					
			### get mRNA sequences ###
			if(missing(mRNA)) {		
				s3utr <- getBM(mart=mart, attributes=c(mRNA_id_type, "3utr"))		
			} else {						
				s3utr <- getBM(mart=mart, attributes=c(mRNA_id_type, "3utr"), 
					filters=mRNA_id_type, values=mRNA)
			}
			
			mRNA <- s3utr[,mRNA_id_type]	
			s3utr <- s3utr$`3utr`	
			names(s3utr) <- mRNA
						
			### get miRNA sequences ###
			if(missing(miRNA)) {
				
				miRNAseq <- getBM(mart=mart, filters=c("biotype"),
					attributes=c(miRNA_id_type, "gene_exon"), values=list("miRNA"))
				
			} else {				
				miRNAseq <- getBM(mart=mart, 
					attributes=c(miRNA_id_type, "gene_exon"),
					filters=c(miRNA_id_type, "biotype"), values=list(miRNA, "miRNA"))
			}
			
			miRNA <- miRNAseq[,miRNA_id_type]	
			miRNAseq <- miRNAseq$gene_exon
			names(miRNAseq) <- miRNA
												
			seedReg <- seedRegions(miRNAseq)
		
			seedReg.rc <- as.character(reverseComplement(DNAStringSet(seedReg)))
		
			mx <- matchSeeds(seedReg.rc, s3utr)
			
			df <- melt(mx)
		
			seedMatrix <- table(df$L2, df$L1)
		}
	}
	
	if(longest3utr) {
		
		geneInfo <- getTranscriptIDwithLongest3UTR(rownames(seedMatrix), mRNA_id_type)
		
		mRNA_id.sel <- geneInfo[, mRNA_id_type]
		
		seedMatrix <- seedMatrix[match(mRNA_id.sel, rownames(seedMatrix)),]
	}
	
	if(returnGeneInfo || convert2genesymbol) {
		
		mRNA_id <- rownames(seedMatrix)
		
		if(!exists("geneInfo")) {				
			
			mart <- useMart(biomart=biomart, dataset=dataset)
			
			geneInfo <- getBM(mart=mart, attributes=c(mRNA_id_type, "ensembl_gene_id", "external_gene_name", "3_utr_start", "3_utr_end"), filters=mRNA_id_type, values=mRNA_id)
			
			geneInfo$utrlen <- abs(geneInfo$`3_utr_end`-geneInfo$`3_utr_start`)
		
			geneInfo <- geneInfo[match(mRNA_id, geneInfo[,mRNA_id_type]),]
		}
		
		if(convert2genesymbol) {
			
			rownames(seedMatrix) <- geneInfo[match(mRNA_id, geneInfo[,mRNA_id_type]),"external_gene_name"]
			
		}		
	} 
	
	if(returnGeneInfo) {
		
		list(seedMatrix=seedMatrix, geneInfo=geneInfo)
		
	} else {
		seedMatrix
	}
}



# Obtain the transcript with longest 3'UTR to represent the genes
getTranscriptIDwithLongest3UTR <- function(mRNA_id, 
	mRNA_id_type="ensembl_transcript_id",
	biomart="ensembl", dataset="hsapiens_gene_ensembl", ...) {			
			
	mart <- useMart(biomart=biomart, dataset=dataset)
	
	geneInfo <- getBM(mart=mart, attributes=c(mRNA_id_type, "ensembl_gene_id", "external_gene_name_id", "3_utr_start", "3_utr_end"), filters=mRNA_id_type, values=mRNA_id)
	
	geneInfo <- geneInfo[!is.na(geneInfo$`3_utr_start`),]
	
	geneInfo$utrlen <- abs(geneInfo$`3_utr_end`-geneInfo$`3_utr_start`)
		
	geneInfo.sel <- geneInfo[match(mRNA_id, geneInfo[,mRNA_id_type]),]
	
	y <- aggregate(data.frame(utrlen=geneInfo.sel$utrlen),
		list(ensembl_gene_id=geneInfo.sel$ensembl_gene_id), max)
	
	y$key <- sprintf("%s:%s", y$ensembl_gene_id, y$utrlen)
	
	geneInfo.sel$key <- sprintf("%s:%s", geneInfo.sel$ensembl_gene_id, geneInfo.sel$utrlen)
		
	geneInfo.sel <- geneInfo.sel[match(y$key, geneInfo.sel$key),]
	
	geneInfo.sel
}

