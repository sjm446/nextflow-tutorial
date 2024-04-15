/* 
 * pipeline input parameters 
 */
params.reads = "$baseDir/data/ggal/gut_{1,2}.fq"
params.transcript = "$baseDir/data/ggal/transcriptome.fa"
params.multiqc = "$baseDir/multiqc"
params.outdir = "results"

log.info """\
         Modified version to create tags
         ===================================
         R N A S E Q - N F   P I P E L I N E    
         ===================================
         transcriptome: ${params.transcript}
         reads        : ${params.reads}
         outdir       : ${params.outdir}
         """
         .stripIndent()

 
/* 
 * define the `index` process that create a binary index 
 * given the transcriptome file
 */
process index {
    resourceLabels region: 'Dublin', user: 'Steve', project: 'project7a'

    input:
    path transcriptome from params.transcript
     
    output:
    path 'index' into index_ch

    script:       
    """
    salmon index --threads $task.cpus -t $transcriptome -i index
    """
}


Channel 
    .fromFilePairs( params.reads, checkIfExists:true )
    .into { read_pairs_ch; read_pairs2_ch } 

/*
 * Run Salmon to perform the quantification of expression using
 * the index and the matched read files
 */
process quantification {
    resourceLabels region: 'Dublin', user: 'Steve', project: 'project7b'

    input:
    path index from index_ch
    tuple val(pair_id), path(reads) from read_pairs_ch
 
    output:
    path(pair_id) into quant_ch
 
    script:
    """
    salmon quant --threads $task.cpus --libType=U -i $index -1 ${reads[0]} -2 ${reads[1]} -o $pair_id
    """
}

/*
 * Run fastQC to check quality of reads files
 */
process fastqc {
    resourceLabels region: 'Dublin', user: 'Steve', project: 'project7c'

    tag "FASTQC on $sample_id"

    input:
    tuple val(sample_id), path(reads) from read_pairs2_ch

    output:
    path("fastqc_${sample_id}_logs") into fastqc_ch

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """  
}

/*
 * Create a report using multiQC for the quantification
 * and fastqc processes
 */
process multiqc {
    resourceLabels region: 'Dublin', user: 'Steve', project: 'project7d'
    publishDir params.outdir, mode:'copy'

    input:
    path('*') from quant_ch.mix(fastqc_ch).collect()

    output:
    path('multiqc_report.html')

    script:
    """
    multiqc .
    """
}

workflow.onComplete { 
	log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong" )
}
