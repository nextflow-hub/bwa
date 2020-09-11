# Nextflow wrapper for `bwa` process.

## Pre-requisites

- Nextflow
- Docker 

**NOTE** If you plan to setup a basic server, then you can refer [minimal-nextflow-setup](https://github.com/nextflow-hub/minimal-nextflow-setup)

## Usage

```
nextflow run https://github.com/nextflow-hub/bwa
```

## Options


- `filePattern`

By default, the process assumes the files to follow the `*_{R1,R2}.fastq.gz` pattern, which could be customized using this option

```
nextflow run https://github.com/nextflow-hub/bwa --filePattern './*_{1,2}.fastq.gz'
```

- `refFasta`

By default, the process assumes the `reference fasta` to be `NC000962_3.fasta`, which could be customized using this option

```
nextflow run https://github.com/nextflow-hub/bwa --refFasta ./XYZ.fasta
```


- `index`

This option can be used to run `bwa index` like so 

```
nextflow run https://github.com/nextflow-hub/bwa --index
```


- `mem`

This option can be used to run `bwa mem` like so 

```
nextflow run https://github.com/nextflow-hub/bwa --mem
```


- `bwaIndexResultsDir`

By default, it stores the result files locally inside the `./results/bwa/index` directory.

```
nextflow run https://github.com/nextflow-hub/bwa --bwaIndexResultsDir ./path/to/custom/bwaIndexResultsDir
```


- `bwaMemResultsDir`

By default, it stores the result files locally inside the `./results/bwa/mem` directory.

```
nextflow run https://github.com/nextflow-hub/bwa --resultsDir ./path/to/custom/bwaMemResultsDir
```


- `samtoolsFaidxResultsDir`

By default, it stores the result files locally inside the `./results/samtools/faidx` directory.

```
nextflow run https://github.com/nextflow-hub/bwa --samtoolsFaidxResultsDir /path/to/custom/samtoolsFaidxResultsDir
```

- `saveMode`

By default, the pipeline publishes the results in the `resultsDir` by copying the relevant output.

You can update this behavior by simply specifying the alternative such as `move` or `link` etc. 

```
nextflow run https://github.com/nextflow-hub/bwa --saveMode move
```

For more information please refer [Nextflow documentation](https://www.nextflow.io/docs/latest/process.html#publishdir)

## Customizing the script

The sole purpose of process wrappers in `nextflow-hub` is to keep the code small, clean and hackable with some basic knowledge of `nextflow` scripting.

If you have specific requirements, you are encouraged to fork/clone and update your version to accomodate your needs. 


## Contribution

Contribution, in all forms, is most welcome!
