#!/usr/bin/env nextflow
OUTDIR = params.outdir+'/'+params.subdir 


csv = file(params.csv)
// Print commit-version of active deployment
//file(params.git)
//    .readLines()
//   .each { println "git commit-hash: "+it }
// Print active container
//container = file(params.container).toRealPath()
//println("container: "+container)


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
    //preprocessing for FastQ files
    

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

/*
process bowtie2_filterfq{
    cpus 16

    when:
        params.filterFastqs
    input:
        set val(sample_id), file(r1), file(r2) from  fq_trimToFilter
    output:
    
        set  val(sample_id), file(r1), file(r2) into fastq_filtered
    script:
        """
        bowtie2 -p ${task.cpus} -q --fr -k 1 --phred33 -t --local --un-conc-gz masked/R%.fastq.gz  -x ${params.RMidx} -1 ${r1} -2 ${r2} -S /dev/null
        """
}
*/

process  create_refidx_hisat2 {
    publishDir "${params.refpath}/hisat2_index", mode :'copy'
    cpus 16
    memory '160 GB'
    //input:
        
    output:
    file "*" into hisat2_ref
   
    
    script:
    """
    hisat2_extract_splice_sites.py $params.gtf > genome.ss 
    hisat2_extract_exons.py $params.gtf > genome.exon
    hisat2_extract_snps_haplotypes_UCSC.py ${params.genome} ${params.snp150} genome
    hisat2-build -p ${task.cpus} --snp genome.snp --haplotype genome.haplotype --exon genome.exon --ss genome.ss ${params.genome} genom_snp_tran
    
    """
}

/*

process  hisat2_align{
/*
--dta/--downstream-transcriptome-assembly
Report alignments tailored for transcript assemblers including StringTie. 
With this option, HISAT2 requires longer anchor lengths for de novo discovery of splice sites. This leads to fewer alignments with short-anchors,
 which helps transcript assemblers improve significantly in computation and memory usage.
--dta-cufflinks
Report alignments tailored specifically for Cufflinks. In addition to what HISAT2 does with the above option (–dta), With this option, HISAT2 looks for novel 
splice sites with three signals (GT/AG, GC/AG, AT/AC), 
but all user-provided splice sites are used irrespective of their signals. HISAT2 produces an optional field, XS:A:[+-], for every spliced alignment.
  
cpus  16
    input:

         set val(sample_id), file(r1), file(r2) from  fastq_trimmed
    output:
        set val(sample_id), file("${sample_id}.trimmed.R1.fq.gz"), file("${sample_id}.trimmed.R2.fq.gz") into fastq_trimmed

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
    -x ${Tidx} \\
    -1 masked/R1.fastq.gz \\
    -2 masked/R2.fastq.gz \\
    -S aligned/alignment.sam
    """
}

process samtools_sort{
    cpus  16
    input:
    output:
    script:
        """
        samtools sort -@ ${task.cpus} -o aligned/alignment.bam -O bam -T aligned/tmp aligned/alignment.sam
        """
}


process  mark_duplicates{ 

    script:
    """
     java -Dpicard.useLegacyParser=false -Xmx${PicardMemory} -jar ${PicardDir}/picard.jar MarkDuplicates \\
         -INPUT aligned/alignment.bam \\
         -OUTPUT aligned/alignment.bam.tmp_picard \\
         -METRICS_FILE aligned/alignment_picardmetrics.csv \\
         -REMOVE_DUPLICATES false \\
         -ASSUME_SORTED true \\
         -MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 2000 \\
         -QUIET true \\
         -VERBOSITY WARNING

    mv -f aligned/alignment.bam.tmp_picard aligned/alignment.bam
    """
}


process samtools_index{

    scripts:
    """
    samtools index -b aligned/alignment.bam aligned/alignment.bai
    """
}

process  stringtie{
    //Expression estimation using protein coding transcripts from GENCODE release 27 as transcriptome model. Novel transcripts are discarded.

    script:
    """
    stringtie -p ${NumThreads} -G ${GTF} -o stringtie/transcript.gtf -A stringtie/gene.tsv -C stringtie/transcript_covered.gtf --rf -B -e bam/alignment.bam
    """


}


*/
