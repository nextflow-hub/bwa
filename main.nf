#!/usr/bin/env nextflow
import java.nio.file.Paths

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

params.mem = false
params.index = false

params.bwaIndexResultsDir = './results/bwa/index'
params.bwaMemResultsDir = './results/bwa/mem'
params.samtoolsFaidxResultsDir = './results/samtools/faidx'

params.saveMode = 'copy'

params.refFasta = "./NC000962_3.fasta"
params.filePattern = "./*_{R1,R2}.fastq.gz"

Channel.value("$workflow.launchDir/$params.refFasta")
        .set { ch_refFasta }

Channel.fromFilePairs(params.filePattern)
        .set { ch_in_bwa }

/*
#==============================================
bwa
#==============================================
*/

process index {
    publishDir params.bwaIndexResultsDir, mode: params.saveMode
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


process mem {
    publishDir params.bwaMemResultsDir, mode: params.saveMode
    container 'quay.io/biocontainers/bwa:0.7.17--hed695b0_7'

    when:
    params.mem

    input:
    path("""${params.bwaIndexResultsDir}""") from Channel.value(Paths.get(params.bwaIndexResultsDir))
    path("""${params.samtoolsFaidxResultsDir}""") from Channel.value(Paths.get(params.samtoolsFaidxResultsDir))
    path refFasta from ch_refFasta
    set genomeFileName, file(genomeReads) from ch_in_bwa

    output:
    file('*.sam') into ch_out_bwaMem


    script:
    // TAG="@RG\\tID:$genomeFileName\\tSM:$genomeFileName\\tLB:$genomeFileName"

    // This alternative TAG hard-codes the PL and PM values. See https://gencore.bio.nyu.edu/variant-calling-pipeline-gatk4/ 
    TAG="@RG\\tID:$genomeFileName\\tLB:$genomeFileName\\tPL:ILLUMINA\\tPM:HISEQ\\tSM:$genomeFileName" 

    """
    cp ${params.bwaIndexResultsDir}/* .
    cp ${params.samtoolsFaidxResultsDir}/* .
    bwa mem -K 100000000 -Y  -R "${TAG}\" ${params.refFasta} ${genomeReads[0]} ${genomeReads[1]} > ${genomeFileName}.sam
    """
}


/*
#==============================================
# extra
#==============================================
*/
