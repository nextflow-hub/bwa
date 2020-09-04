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

Channel.fromFilePairs(params.filePattern)
        .into { ch_in_bwa }

Channel.value("$workflow.launchDir/$params.refFasta")
        .set { ch_refFasta }

/*
#==============================================
bwa
#==============================================
*/

process bwa {
    publishDir params.resultsDir, mode: params.saveMode
    container 'bwa'


    input:
    set genomeFileName, file(genomeReads) from ch_in_bwa

    output:
    path bwa into ch_out_bwa


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
