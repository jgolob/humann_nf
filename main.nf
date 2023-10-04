#!/usr/bin/env nextflow

/*
  humann_nf: A pipeline to wrap around and run Humann3.

*/

// Using DSL-2
nextflow.enable.dsl=2

// containers
container__humann = 'quay.io/biocontainers/humann:3.8--pyh7cba7a3_0'
container__vsearch = 'quay.io/biocontainers/vsearch:2.17.0--h95f258a_1'

// Default parameters
params.manifest = false
params.help = false
params.chocophlan = false
params.uniref = false
params.metaphlan = false
params.metaphlan_index = 'mpa_vOct22_CHOCOPhlAnSGB_202212'

params.output = './output'

workflow Humann3 {
    take:
        paired_ch
        se_ch
        chocophlan
        uniref
        metaphlan
        metaphlan_index

    main:

    // merge the PE reads
    Merge_pairs(
        paired_ch.map{
            [it['specimen'], file(it['R1']), file(it['R2'])]
        }
    )

    // Run Humann3
    Run_Humann(
        se_ch.mix(
            Merge_pairs.out,
        ),
        chocophlan,
        uniref,
        metaphlan,
        metaphlan_index
    )

}

process Merge_pairs {
    container "${container__vsearch}"
    label 'io_limited'
    errorStrategy 'ignore'

    input:
    tuple val(specimen), path(R1), path(R2)

    output:
    tuple val(specimen), path("${specimen}.fastq.gz")

    """
    vsearch --fastq_mergepairs \
    ${R1} --reverse ${R2} \
    --fastqout ${specimen}.fastq

    gzip ${specimen}.fastq
    """
}

process Run_Humann {
    container "${container__humann}"
    label 'multithread'
    errorStrategy 'ignore'
    publishDir "${params.output}", mode: 'copy'

    beforeScript 'ulimit -Ss unlimited'

    input:
    tuple val(specimen), path(R1)
    path(chocophlan), stageAs: 'refs/chocophlan/'
    path(uniref), stageAs: 'refs/uniref/'
    path(metaphlan), stageAs: 'refs/metaphlan/'
    val(metaphlan_index)

    output:
    tuple val(specimen), path("out/${specimen}_genefamilies.tsv"), path("out/${specimen}_pathabundance.tsv"), path("out/${specimen}_pathcoverage.tsv")
    """
    humann3 --input ${R1} \
    --protein-database ${uniref} \
    --nucleotide-database ${chocophlan} \
    --metaphlan-options "-t rel_ab  --index ${metaphlan_index} --bowtie2db ${metaphlan} --offline --force" \
    --threads ${task.cpus} \
    --output out/
    """
}



def helpMessage() {
    log.info"""
    Workflow to run Humann3 on metagenomics specimens.

    ** Reads must be preprocessed first: filtered, trimmed, and with human reads removed.
    Usage:
    --manifest      Path to a read-pairs to be analyzed (REQUIRED)
                    'specimen'                  This can be any string, which is a sequence of alpha-numeric characters, underscores, 
                                                or dashes and no spaces, that is less than 64 characters.
                                                MUST be unique for this data set. 
                    'R1'                        First read, in fasta format and gzipped. Must be preprocessed.
                    'R2'                        Second / reverse read, in fasta format and gzipped  (optional)
    --chocophlan    Path to chocophlan references directory (REQUIRED)
    --uniref        Path to uncompressed uniref diamond reference directory (REQUIRED)
    --metaphlan     Path to metaphlan reference directory (REQUIRED)
    --metaphlan_index   Current index (default: mpa_vOct22_CHOCOPhlAnSGB_202212)
    --output        Path where to place the output files
    
    Parameters:
    """.stripIndent()
}

def ReadManifest(manifest_file){
    manifest_file.splitCsv(
        header: true, 
        sep: ","
    ).branch{
        valid_paired:  (it.specimen != null) && (it.R1 != null ) && (it.R1 != "" ) && (!file(it.R1).isEmpty()) && (it.R2 != null ) && (it.R2 != "") && (!file(it.R2).isEmpty())
        valid_single:  (it.specimen != null) && (it.R1 != null ) && (it.R1 != "" ) && (!file(it.R1).isEmpty())
        other: true
    }
}

workflow {
    if ((!params.chocophlan) | (!params.uniref) | (!params.manifest) | (!params.metaphlan) | params.help) {
        helpMessage()
        exit 0
    }

    // Load manifest!
    manifest = ReadManifest(
        Channel.from(
            file(params.manifest)
        )
    )

    Humann3(
        manifest.valid_paired,
        manifest.valid_single,
        params.chocophlan,
        params.uniref,
        params.metaphlan,
        params.metaphlan_index
    )

}
