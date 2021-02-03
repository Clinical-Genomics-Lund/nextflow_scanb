OUTDIR = params.outdir+'/'+params.subdir
process star_alignment{
    tag "${smpl_id}"
    container = "/fs1/resources/containers/rnaseqfus_active.sif"
    errorStrategy 'ignore'
    publishDir "$OUTDIR/bam", mode :'copy'
    cpus = 16
    memory 40.GB
    //when:
    //        params.qc || params.star

    input:
            tuple val(smpl_id) , file(read1), file(read2) 

    output:
            tuple val(smpl_id), file("${smpl_id}.Aligned.sortedByCoord.out.bam"), emit: star_bam
            tuple val(smpl_id), file("${smpl_id}.Log.final.out"), emit: star_log        
    script: 
    """
    STAR --genomeDir ${params.ref_genome_dir} \\
            --readFilesIn ${read1} ${read2} \\
            --runThreadN ${task.cpus} \\
            --outSAMtype BAM SortedByCoordinate \\
            --readFilesCommand zcat \\
            --limitBAMsortRAM 10000000000 

    mv Aligned.sortedByCoord.out.bam  ${smpl_id}.Aligned.sortedByCoord.out.bam
    mv Log.final.out ${smpl_id}.Log.final.out
    """
}

process index_bam {
        tag "${smpl_id}"
        container = "/fs1/resources/containers/rnaseqfus_active.sif"
        publishDir "$OUTDIR/bam", mode :'copy'

        input:
                tuple val(smpl_id), file(bam)
        output:
                tuple val(smpl_id), file("${smpl_id}.Aligned.sortedByCoord.out.bam"), file("${smpl_id}.Aligned.sortedByCoord.out.bam.bai") , emit: bam_bai
                
        script:
        """
        sambamba index --show-progress -t 8 ${bam}
        """     
}  

process geneBodyCov{
    tag "${smpl_id}"
	publishDir "$OUTDIR/qc", mode:'copy'
	//errorStrategy 'ignore'
	cpus 1
    container = "/fs1/resources/containers/rnaseqfus_active.sif"
	
	//when:
	//	params.qc || params.bodyCov
	
	input:
	
		tuple val(smpl_id),file(bam_f)
				
	output:
		
		tuple val(smpl_id), file("${smpl_id}.geneBodyCoverage.txt"), emit: geneBodyCov
		
	script:
	
	"""
	samtools view  -s 0.3 -b ${bam_f} -o ${smpl_id}.subsample.bam
	sambamba index --show-progress -t 8 ${smpl_id}.subsample.bam
	geneBody_coverage.py -i ${smpl_id}.subsample.bam -r ${params.ref_rseqc_bed} -o ${smpl_id}
	"""
	}

process provider{
    container = "/fs1/resources/containers/rnaseqfus_active.sif"
	tag "${smpl_id}"
	//publishDir "$OUTDIR/qc" , mode:'copy'
	//errorStrategy 'ignore'

	input:
		tuple val(smpl_id), file(bam_f) 
		
	
	output:
		tuple val(smpl_id), file("${smpl_id}.genotypes"), emit: provider_geno

	script:

	prefix = "${smpl_id}"
	"""
	provider.pl  --out ${prefix} --bed ${params.ref_bed} --bam ${bam_f} --bedxy ${params.ref_bedXy}
	"""
	}
	
process postaln_qc_rna {
	//publishDir "$OUTDIR/finalResults" , mode:'copy'
    container = "/fs1/resources/containers/rnaseqfus_active.sif"
	tag "${smpl_id}"

	input:
		tuple val(smpl_id), file(star_log), file(geneCov), file(provider_geno)
	
	output:
		tuple val(smpl_id), file("${smpl_id}.STAR.rnaseq_QC"), emit: rnaseq_qc 

	script:
	"""
	postaln_qc_rna.R \\
		--star ${star_log} \\
		--id '${smpl_id}' \\
		--provider ${provider_geno} \\
		--genebody ${geneCov}> '${smpl_id}.STAR.rnaseq_QC'
	"""
	} 

process  register_to_cdm{
	publishDir "${params.crondir}/qc", mode: 'copy', overwrite: true
	cpus 1
	memory '8 GB'
	time '1h'
	input:
		tuple val(smpl_id), file(qc), fastq_r1, fatq_r2 ,clarity_id, pool_id  
	 output:
		tuple val(smpl_id), file("${smpl_id}.cdm")
	 script:
	 
	 parts = fastq_r1.toString().split('/')
	 parts.println()
	 idx= parts.findIndexOf {it ==~ /......_......_...._........../}
	 rundir= parts[0..idx].join("/")
	 
	 
	"""
	echo "--run-folder ${rundir} --sample-id ${smpl_id} --assay rnaseq-fusion --qc $OUTDIR/finalResults/${qc}" > ${smpl_id}.cdm 
	"""

}
