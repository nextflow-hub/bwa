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

process bwa {
    publishDir params.resultsDir, mode: params.saveMode
    container 'quay.io/biocontainers/bwa:0.7.17--hed695b0_7'


    input:
    set genomeFileName, file(genomeReads) from ch_in_bwa

    output:
    tuple file('*.amb'),
            file('*.ann'),
            file('*.bwt'),
            file('*.fai'),
            file('*.pac'),
            file('*.sa') into ch_out_bwa


    script:
    genomeName = genomeFileName.toString().split("\\_")[0]

    """
    bwa index $params.refFasta
    """
}


/*
#==============================================
# extra
#==============================================
*/
