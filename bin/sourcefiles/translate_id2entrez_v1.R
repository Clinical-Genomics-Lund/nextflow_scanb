
###########################################
#	Johan Vallon-Christersson             #
#	johan.vallon-christersson@med.lu.se   #
###########################################


####
# function: translate_id2entrez
#
#	INPUT:
#		- vector of gene id to translate (id of type used in our StringTie pipeline, i.e., 'Gene.ID')
#		- gene annotation object uncluding columns 'Gene.ID' and 'EntrezGene' (the gene annotation file created for our StringTie pipeline).
#		- vector of unique gene id to return data for (genes to be included in returned matrix). Must be same type of id found in column 'Gene.ID' in gene.tsv file.
#	OUTPUT:
#		- vector of EntrezGene or vector with EntrezGeneMust with appended prefix 'e' (as required by our SSP models).
#
#	v1:	first implementation
#		Returned vector will have NA for id not found in gene annotation object.
#		Note that in gene annotation object all id (Gene.ID) without an EntrezGene have character "NA" in column EntrezGene and these will be returned
# 		



		#	# examples and manual stuff
		#				
		#		# load gene annotation
		#			load("Gene.ID.ann.Rdata")
		#				
		#			# some gene id to translate
		#			some.gene.id <-  Gene.ID.ann$Gene.ID[c(1:5)]
		#				# some.gene.id <-  Gene.ID.ann$Gene.ID[c(1:5, 19625)]
		#				# some.gene.id <- c(some.gene.id, "gaga", some.gene.id)
		#			
		#		# run function
		#			myid <- translate_id2entrez(id=some.gene.id, ann=Gene.ID.ann, e=TRUE, report=TRUE)
		#
		#				head(myid)
			

#####################################
# function translate_id2entrez
#####################################

	# function
	translate_id2entrez <- function(id, ann, e=FALSE, report=FALSE){
					
		# find id index in sum_fpkm_gene$Gene.ID
			id.i <- match(id, ann$Gene.ID)
		
		if(report){
			# test verify that all provided id are found
			if(length(which(is.na(id.i)==TRUE))==0){
				cat("all id found in ann, ")				
			}else{
				cat(length(which(is.na(id.i)==TRUE)), "id not found in ann, ")
			}
		}

		# get EntrezGene and assign id as names 
			entrez <- ann$EntrezGene[id.i]
			names(entrez) <- id
			
			entrezNA <- length(which(entrez[is.na(entrez)==FALSE]=="NA"))
			
			if(report){
				if(entrezNA>0){
					cat(entrezNA, "found id have EntrezGene NA, ")
				}else{
					cat("all found id have EntrezGene, ")
				}				
			}

		# add prefix	
			if(e){
				if(report){
					cat("adding prefix e, ")				
				}
				entrez[is.na(entrez)==FALSE] <- paste("e", entrez[is.na(entrez)==FALSE], sep="")
			}				
			
		if(report){
			cat("done")	
		}

	return(entrez)
		
	}# end function
	
#####################################
#####################################
