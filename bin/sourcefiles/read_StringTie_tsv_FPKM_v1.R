
###########################################
#	Johan Vallon-Christersson             #
#	johan.vallon-christersson@med.lu.se   #
###########################################


####
# function: read_StringTie_tsv_FPKM
#
#	INPUT:
#		- gene.tsv file from StingTie. Must contain column holding gene id with column name 'Gene ID' or 'Gene.ID' and column 'FPKM' with geneexpression data.
#		- vector of unique gene id to return data for (genes to be included in returned matrix). Must be same type of id found in column 'Gene.ID' in gene.tsv file.
#	OUTPUT:
#		- returns genematrix with one column (FPKM data sum on gene id) and rows for gene id (row names are id)

#	v1:	first implementation
# 		the gene.tsv is read using read.delim() and invalid characters in column names are translated to ".". That is, column name 'Gene ID' will be translated to 'Gene.ID'.



		#	# examples and manual stuff
		#			
		#		# specify input gene.tsv 
		#			gene.tsv <- "gene.tsv"
		#				
		#		# load gene annotation
		#			load("Gene.ID.ann.Rdata")
		#				
		#			# some gene id to collect
		#			some.gene.id <-  Gene.ID.ann$Gene.ID[1:5]
		#			
		#		# run function
		#			mymatrix <- read_StringTie_tsv_FPKM(tsv=gene.tsv, id=some.gene.id, report=TRUE)
		#
		#				head(mymatrix)
			

#####################################
# function read_StringTie_tsv_FPKM
#####################################

	# function
	read_StringTie_tsv_FPKM <- function(tsv, id, report=FALSE){
				
			# number of columns and rows in genematrix 
				my.ncol <- 1
				my.nrow <- length(id)
			
				# test verify that id are unique  
				if(length(id)==length(unique(id))){
					if(report){
						cat("provided id(s) are unique, ")	
					}
				}else{
					stop("mfkr id")
				}			
			
			# create genematrix
			genematrix <- matrix(nrow=my.nrow, ncol=my.ncol)
			
				# rownames and colname
					row.names(genematrix) <- id
					colnames(genematrix) <- "FPKMsum"
				
					# head(genematrix)

				# read gene.tsv
					if(report){
						cat("creating matrix, ")
					}
						
					# read from file
						read.data <- read.delim(file=tsv, as.is=TRUE)
							# head(read.data)
						
					# aggregate fpkm on gene 	
						sum_fpkm_gene <- aggregate(FPKM ~ Gene.ID, data=read.data, sum)
							# head(sum_fpkm_gene)

					# find id index in sum_fpkm_gene$Gene.ID
						id.i <- match(id, sum_fpkm_gene$Gene.ID)

					# test verify that all provided id are found
						if(length(which(is.na(id.i)==TRUE))==0){
							# add id.i FPKM to genematrix
								genematrix[,1] <- sum_fpkm_gene$FPKM[id.i]
						}else{
							stop("provided id not found in tsv")
						}
	
					if(report){
						cat("done")	
					}
				return(genematrix)
		
	}# end function
	
#####################################
#####################################
