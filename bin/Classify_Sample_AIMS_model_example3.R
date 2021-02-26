
###########################################
#	Johan Vallon-Christersson             #
#	johan.vallon-christersson@med.lu.se   #
###########################################


#####################################
# example of using function applySSP()
#####################################

	# the eaxamples below is for running in the parent directory of directories 'testdata' and 'sourcefiles'

	# example applySSP_v1.2 with txt and plot outpur


	# 1 #	specify infiles and source direcory
					
		# specify your input gene.tsv to classify
			tsv <- "testdata/gene.tsv"
		
		# specify your input ssp model to use for classification 
			ssp <- "testdata/Training_Run19081Genes_noNorm_SSP.subtypeMost.Fcc15_5x5foldCV.num.rules.50_24.selRules.AIMS.GS.RData"
				
		# specify sourcefiles folder containing a number of required files including sourcefile for function applySSP()
			source <- "sourcefiles"


	# 2 #	 example script for running a classification

		# load the  applySSP() funtion 
			source(paste(source, "applySSP_v1.2.R", sep="/"))
										
		# run function
			myresults <- applySSP(tsv, ssp, source, plot=FALSE, txt=TRUE, add.is.num=TRUE)
		
			print(myresults)
			


	# 3 # test with other SSP

		ssp <- "testdata/Training_Run19081Genes_noNorm_SSP.PAM50subtype4Most.Fcc15_5x5foldCV.num.rules.50_21.selRules.AIMS.GS.RData"
		
		# run function
			myresults <- applySSP(tsv, ssp, source)
		
			print(myresults)


	# 4 # test with yet another SSP

		ssp <- "testdata/Training_Run19081Genes_noNorm_SSP.ER_v2.Fcc15_5x5foldCV.num.rules.50_19.selRules.AIMS.GS.RData"
		
		# run function
			myresults <- applySSP(tsv, ssp, source)
		
			print(myresults)



	# 5 # test with yet another SSP

		ssp <- "testdata/Training_Run19081Genes_noNorm_SSP.scaled.ROR.tot.asT0.c005.Fcc15_5x5foldCV.num.rules.50_21.selRules.AIMS.GS.RData"
		
		# run function
			myresults <- applySSP(tsv, ssp, source)
		
			print(myresults)
