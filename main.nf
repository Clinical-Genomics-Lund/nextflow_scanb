#!/usr/bin/env nextflow

OUTDIR = params.outdir+'/'+params.subdir 

csv = file(params.csv)
//Print commit-version of active deployment

file(params.git)
    .readLines()
   .each { println "git commit-hash: "+it }
// Print active container
container = file(params.container).toRealPath()
println("container: "+container)


nextflow.preview.dsl = 2
include  './modules/aln_qc.nf' params(params)
include  './modules/expr_scanb.nf' params(params)

workflow.onComplete {

	def msg = """\
		Pipeline execution summary
		---------------------------
		Completed at: ${workflow.complete}
		Duration    : ${workflow.duration}
		Success     : ${workflow.success}
		scriptFile  : ${workflow.scriptFile}
		workDir     : ${workflow.workDir}
		exit status : ${workflow.exitStatus}
		errorMessage: ${workflow.errorMessage}
		errorReport :
		"""
		.stripIndent()
	def error = """\
		${workflow.errorReport}
		"""
		.stripIndent()

	base = csv.getBaseName()
	logFile = file("/fs1/results/cron/logs/" + base + ".complete")
	logFile.text = msg
	logFile.append(error)
}
Channel
    .fromPath(params.csv)
    .splitCsv(header:true)
    .map{row -> tuple(row.id,row.read1, row.read2,row.clarity_sample_id, row.clarity_pool_id)}
    .set{cdm_meta}	

Channel
    .fromPath(params.csv)
    .splitCsv(header:true)
    .map{ row-> tuple(row.id, file(row.read1), file(row.read2)) }
    .set {reads}

Channel
    .fromPath(params.csv)
    .splitCsv(header:true)
    .map{ row-> tuple(row.id, file(row.read1), file(row.read2)) }
    .set {fastqs}

workflow qc_workflow{
	take: fastqs
	main:
		star_alignment(fastqs)
		index_bam(star_alignment.out.star_bam)
		geneBodyCov(star_alignment.out.star_bam)
		provider(star_alignment.out.star_bam)
		qc_files = star_alignment.out.star_log.join(geneBodyCov.out.geneBodyCov.join(provider.out.provider_geno))
		postaln_qc_rna(qc_files)
		cdm_data=postaln_qc_rna.out.rnaseq_qc.join(cdm_meta)
		cdm_data.view()
		register_to_cdm(cdm_data)
	emit:
	register_to_cdm.out
	
}


workflow expr_workflow{
	 take:fastqs
	 main:
		fastp(fastqs)
		bowtie2_filterfq(fastp.out.fq_trimToFilter)
		//create_refidx_hisat2()
		hisat2_align(bowtie2_filterfq.out.fastq_filtered)
		samtools_sort(hisat2_align.out.hs2_sam)
		mark_duplicates(samtools_sort.out.bam)
		samtools_index(mark_duplicates.out.bam_dup)
		stringtie(samtools_index.out.bambai)
	emit:
	stringtie.out
}

workflow{
	take:fastqs
	main:
	fastqs.view()
	expr_workflow()
	qc_workflow()
    
	}
	
