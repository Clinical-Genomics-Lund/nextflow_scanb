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
    .map{ row-> tuple(row.id, file(row.read1), file(row.read2)) }
    .set {reads}


process  fastp{
	tag "${sample_id}" 
    //preprocessing for FastQ files
    publishDir "$OUTDIR/qc", mode :'copy'
    
    input:
        set val(sample_id), file(r1), file(r2) from reads

    output:
        set val(sample_id), file("${sample_id}.trimmed.R1.fq.gz"), file("${sample_id}.trimmed.R2.fq.gz") into fastq_trimmed, fq_trimToFilter
        set val(sample_id), file("${sample_id}.fastp.json"), file("${sample_id}.fastp.html") into fastp_report
    
    script:
    """
    fastp -i $r1 -I $r2 -o  ${sample_id}.trimmed.R1.fq.gz -O ${sample_id}.trimmed.R2.fq.gz -j fastp.json -h fastp.html
    mv  fastp.json ${sample_id}.fastp.json
    mv  fastp.html  ${sample_id}.fastp.html 
    """

}

process bowtie2_filterfq {
	tag "${sample_id}" 
    cpus 16
    memory 40.GB 
    
    when:
        params.filterFastqs
    input:
        set val(sample_id), file(r1), file(r2) from  fq_trimToFilter
    output:
    
        set val(sample_id), file("${sample_id}.trimmed.filtered.R1.fastq.gz"), file("${sample_id}.trimmed.filtered.R2.fastq.gz") into fastq_filtered

    script:
        """
        bowtie2 -p ${task.cpus} -q --fr -k 1 --phred33 -t --local --un-conc-gz ${sample_id}.trimmed.filtered.R%.fastq.gz -x ${params.RMidx} -1 ${r1} -2 ${r2} -S /dev/null
        """
}


process create_refidx_hisat2{
	tag "${sample_id}" 
    publishDir "${params.refpath}/hisat2_index", mode :'copy'
    cpus 16
    memory '160 GB'
    
    output:
        file "*" into hisat2_ref
    
    when:
        params.ht2ref
    
    script:
    """
    hisat2_extract_splice_sites.py $params.gtf > genome.ss 
    hisat2_extract_exons.py $params.gtf > genome.exon
    hisat2_extract_snps_haplotypes_UCSC.py ${params.genome} ${params.snp150} genome
    hisat2-build -p ${task.cpus} --snp genome.snp --haplotype genome.haplotype --exon genome.exon --ss genome.ss ${params.genome} genome_snp_tran
    """
}

process  hisat2_align {
	tag "${sample_id}" 
    cpus 16
    memory '64 GB'

    input:
        set val(sample_id), file(r1), file(r2) from  fastq_filtered
    
    output:
        set val(sample_id), file("${sample_id}.hst2Aligned.sam") into hs2_sam

    script:
    """
    hisat2 \\
        -p ${task.cpus} \\
        -q --fr --phred33 -t --dta --dta-cufflink --new-summary --no-unal --non-deterministic \\
        --novel-splicesite-outfile aligned/splicesites.tsv \\
        --rna-strandness RF \\
        --summary-file aligned/summary.txt \\
        --rg PL:Illumina \\
        --rg CN:SCANB-prim \\
        --rg-id SAMPLE.l.r.n.m.c.lib.g2 \\
        --rg PU:HMTL5BGXF \\
        --rg SM:SAMPLE \\
        --rg LB:SAMPLE.l.r.n.m.c.lib \\
        -x ${params.indxDir}/genom_snp_tran \\
        -1 ${r1} \\
        -2 ${r2} \\
        -S ${sample_id}.hst2Aligned.sam
    """
}

process samtools_sort{
	tag "${sample_id}" 
    publishDir "$OUTDIR/bam", mode :'copy'
    cpus  16
    
    input:
        set val(sample_id), file(sam) from  hs2_sam
    
    output:
        set val(sample_id), file("${sample_id}.hst2Aligned.sorted.bam") into bam

    script:
    """
    samtools sort -@ ${task.cpus} -o ${sample_id}.hst2Aligned.sorted.bam -O bam -T ./tmp ${sam}
    """
}


process  mark_duplicates{
	tag "${sample_id}" 
    container="/fs1/resources/containers/wgs_2020-03-25.sif"
    publishDir "$OUTDIR/bam", mode :'copy'
    cpus 16
    memory '32 GB'

    
    input:
        set val(sample_id), file(bamfile) from bam
    
    output:
        set val(sample_id), file("${sample_id}.markdup.bam") into bam_dup
    
    script:
    """
    java -Dpicard.useLegacyParser=false -Xms24G -Xmx24G -jar /opt/conda/envs/CMD-WGS/share/picard-2.21.2-1/picard.jar MarkDuplicates \\
         -INPUT ${bamfile} \\
         -OUTPUT ${sample_id}.alignment.bam.tmp_picard \\
         -METRICS_FILE ${sample_id}.alignment_picardmetrics.csv \\
         -REMOVE_DUPLICATES false \\
         -ASSUME_SORTED true \\
         -MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 2000 \\
         -QUIET true \\
         -VERBOSITY WARNING
    mv  -f  ${sample_id}.alignment.bam.tmp_picard  ${sample_id}.markdup.bam
    
    """
}

process samtools_index{
	tag "${sample_id}" 
    publishDir "$OUTDIR/bam", mode :'copy'
    
    input:
        set val(sample_id), file(bam) from bam_dup
    
    output:
        set val(sample_id), file("${sample_id}.markdup.bam"), file("${sample_id}.markdup.bam.bai") into bambai
    
    script:
    """
    samtools index -b $bam ${bam}.bai
    """
}

process  stringtie {
	tag "${sample_id}" 
    cpus 16
    memory '32 GB'
    
    //Expression estimation using protein coding transcripts from GENCODE release 27 as transcriptome model. Novel transcripts are discarded.
    publishDir "$OUTDIR/stringtie", mode :'copy' 
    
    input:
        set val(sample_id), file(bam), file(bai) from  bambai

    output:
        set val(sample_id), file("${sample_id}.out.gtf"), file("${sample_id}.gene_abund.tsv"), file("${sample_id}.transcript_covered.gtf") into abund_ch

    script:
    """
    stringtie -p ${task.cpus} -G ${params.gtf} -o ${sample_id}.out.gtf -A ${sample_id}.gene_abund.tsv -C ${sample_id}.transcript_covered.gtf --rf -B -e ${bam}
    """
}
