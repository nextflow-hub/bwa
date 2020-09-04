#!/usr/bin/env nextflow

/*
#==============================================
code documentation
#==============================================
*/


/*
#==============================================
params
#==============================================
*/

params.resultsDir = 'results/bwa'
params.saveMode = 'copy'
params.filePattern = "./*_{R1,R2}.fastq.gz"
params.mem = false
params.index = false
params.refFasta = "NC000962_3.fasta"

Channel.value("$workflow.launchDir/$params.refFasta")
        .set { ch_refFasta }

Channel.fromFilePairs(params.filePattern)
        .set { ch_in_bwa }

/*
#==============================================
bwa
#==============================================
*/

process bwaIndex {
    publishDir params.resultsDir, mode: params.saveMode
    container 'quay.io/biocontainers/bwa:0.7.17--hed695b0_7'

    when:
    params.index

    input:
    path refFasta from ch_refFasta

    output:
    tuple file('*.amb'),
            file('*.ann'),
            file('*.bwt'),
            file('*.pac'),
            file('*.sa') into ch_out_bwaIndex


    script:

    """
    bwa index $params.refFasta
    """
}


process bwaMem {
    publishDir params.resultsDir, mode: params.saveMode
    container 'quay.io/biocontainers/bwa:0.7.17--hed695b0_7'

    when:
    params.mem

    input:
    path refFasta from ch_refFasta
    set genomeFileName, file(genomeReads) from ch_in_bwa

    output:
    file('*.bam') into ch_out_bwaMem


    script:
    TAG="@RG\\tID:$genomeFileName\\tSM:$genomeFileName\\tLB:$genomeFileName"

    """
    bwa mem -R \"${TAG}\" ${params.refFasta} ${genomeReads[0]} ${genomeReads[1]}
    """
}


/*
#==============================================
# extra
#==============================================
*/
