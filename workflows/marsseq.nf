/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_marsseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MARSSEQ {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    
    main:

    ch_versions             = Channel.empty()
    ch_multiqc_files        = Channel.empty()
    ch_fasta                = file(params.fasta, checkIfExists: true)
    ch_gtf                  = file(params.gtf, checkIfExists: true)
    ch_bowtie_index         = file(params.bowtie2_index, checkIfExists: true)
    ch_star_index           = file(params.star_index, checkIfExists: true)
    ch_ercc_regions         = Channel.fromPath("$projectDir/data/ercc-regions.tsv")
    ch_oligos               = Channel.fromPath("$projectDir/data/oligos.txt")
    ch_spike_seq            = Channel.fromPath("$projectDir/data/spike-seq.txt")
    ch_spike_concentrations = Channel.fromPath("$projectDir/data/spike-concentrations.txt")
    
    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    // PREPARE_PIPELINE (
    //     INPUT_CHECK.out.reads.map { it[0].amp_batches },
    //     INPUT_CHECK.out.reads.map { it[0].seq_batches },
    //     INPUT_CHECK.out.reads.map { it[0].well_cells },
    //     ch_gtf,
    //     ch_ercc_regions,
    //     INPUT_CHECK.out.reads
    // )
    // ch_versions = ch_versions.mix(PREPARE_PIPELINE.out.versions)

    // LABEL_READS (
    //     ch_oligos,
    //     PREPARE_PIPELINE.out.amp_batches,
    //     PREPARE_PIPELINE.out.seq_batches,
    //     PREPARE_PIPELINE.out.reads
    // )
    // ch_versions = ch_versions.mix(LABEL_READS.out.versions)

    // ALIGN_READS ( LABEL_READS.out.read, ch_bowtie_index, LABEL_READS.out.qc )
    // ch_versions = ch_versions.mix(ALIGN_READS.out.versions)

    // // merge sam files into one file
    // ch_aligned_reads = ALIGN_READS.out.reads
    //     .map { meta, sam -> [ meta.id, sam ] }
    //     .groupTuple(by: [0], sort: { it.name })
    //     .map { batch, sams -> [ [ "id": batch ], sams ] }

    // // merged aligned SAM files
    // MERGE_READS ( ch_aligned_reads )
    // ch_versions = ch_versions.mix(MERGE_READS.out.versions)

    // DEMULTIPLEX_READS (
    //     MERGE_READS.out.file_out,
    //     PREPARE_PIPELINE.out.amp_batches,
    //     PREPARE_PIPELINE.out.seq_batches,
    //     PREPARE_PIPELINE.out.wells_cells,
    //     PREPARE_PIPELINE.out.gene_intervals,
    //     ch_spike_seq,
    //     ch_spike_concentrations,
    //     ch_oligos
    // )
    // ch_versions = ch_versions.mix(DEMULTIPLEX_READS.out.versions)

    // QC_REPORT (
    //     DEMULTIPLEX_READS.out.qc_rd.map { meta, rds -> [ ["id": meta.id], rds ] }.groupTuple(),
    //     DEMULTIPLEX_READS.out.qc_pdf.map { meta, pdf -> [ ["id": meta.id], pdf ] }.groupTuple(),
    //     PREPARE_PIPELINE.out.amp_batches,
    //     PREPARE_PIPELINE.out.wells_cells
    // )
    // ch_versions = ch_versions.mix(QC_REPORT.out.versions)

    // //
    // // MODULE: Velocity
    // //
    // if (params.velocity) {
    //     VELOCITY (
    //         INPUT_CHECK.out.reads.map { it[0].amp_batches },
    //         INPUT_CHECK.out.reads.map { it[0].well_cells },
    //         PREPARE_PIPELINE.out.reads,
    //         ch_star_index,
    //         ch_gtf
    //     )
    //     ch_versions = ch_versions.mix(VELOCITY.out.versions)

    //     ch_multiqc_files = ch_multiqc_files.mix(VELOCITY.out.catadapt_multiqc.collect{it[1]}.ifEmpty([]))
    //     ch_multiqc_files = ch_multiqc_files.mix(VELOCITY.out.star_multiqc.collect{it[1]}.ifEmpty([]))
    // }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  + 'pipeline_software_' +  'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions = ch_versions                            // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
