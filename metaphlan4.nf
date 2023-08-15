#!/usr/bin/env nextflow

/*
  A pipeline to just run metaphan (v4) on specimens

*/

// Using DSL-2
nextflow.enable.dsl=2

// containers
container__metaphlan = 'golob/metaphlan:v4.06'
container__vsearch = 'quay.io/biocontainers/vsearch:2.17.0--h95f258a_1'

// Default parameters
params.manifest = false
params.help = false
params.metaphlan = false
params.metaphlan_index = 'mpa_vOct22_CHOCOPhlAnSGB_202212'

params.output = './output'

workflow MetaPhlAn4 {
    take:
        paired_ch
        se_ch
        metaphlan
        metaphlan_index

    main:

    println(
        metaphlan_index
    )
    // merge the PE reads
    Merge_pairs(
        paired_ch.map{
            [it['specimen'], file(it['R1']), file(it['R2'])]
        }
    )

    // Run MetaPhlAn4
    Run_MetaPhlAn4(
        se_ch.mix(
            Merge_pairs.out,
        ),
        metaphlan,
        metaphlan_index
    )

    Run_MetaPhlAn4.out.view()


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

process Run_MetaPhlAn4 {
    container "${container__metaphlan}"
    label 'mid_memory'
    errorStrategy 'ignore'
    publishDir "${params.output}", mode: 'copy'

    beforeScript 'ulimit -Ss unlimited'

    input:
        tuple val(specimen), path(R1)
        path(metaphlan)
        val(metaphlan_index)

    output:
        tuple val (specimen), path("${specimen}.mp4_rel_ab_w_read_stats.tsv"), path("${R1}.bowtie2out.txt")
    
    """
    metaphlan ${R1} ${specimen}.mp4_rel_ab_w_read_stats.tsv \
    --input_type fastq \
    -t rel_ab_w_read_stats \
    --offline \
    --index "${metaphlan_index}" \
    --bowtie2db ${metaphlan}/ \
    --nproc ${task.cpus} \
    --sample_id ${specimen} \
    --sample_id_key ${specimen} \
    --output out/
    """
}



def helpMessage() {
    log.info"""
    Workflow to run MetaPhlAn4 on metagenomics specimens.

    ** Reads must be preprocessed first: filtered, trimmed, and with human reads removed.
    Usage:
    --manifest      Path to a read-pairs to be analyzed (REQUIRED)
                    'specimen'                  This can be any string, which is a sequence of alpha-numeric characters, underscores, 
                                                or dashes and no spaces, that is less than 64 characters.
                                                MUST be unique for this data set. 
                    'R1'                        First read, in fasta format and gzipped. Must be preprocessed.
                    'R2'                        Second / reverse read, in fasta format and gzipped  (optional)
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
    if ((!params.manifest) | (!params.metaphlan) | params.help) {
        helpMessage()
        exit 0
    }

    // Load manifest!
    manifest = ReadManifest(
        Channel.from(
            file(params.manifest)
        )
    )

    MetaPhlAn4(
        manifest.valid_paired,
        manifest.valid_single,
        params.metaphlan,
        params.metaphlan_index
    )

}
