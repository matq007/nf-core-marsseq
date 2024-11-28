//
// Prepare reference genome files
//

include { ERCC_CREATE          } from '../../modules/local/ercc/main'
include { CAT_CAT as CAT_FASTA } from '../../modules/nf-core/cat/cat/main'
include { BOWTIE2_BUILD        } from '../../modules/nf-core/bowtie2/build/main'
include { STAR_GENOMEGENERATE  } from '../../modules/nf-core/star/genomegenerate/main'

workflow PREPARE_GENOME {

    take:
    aligner         // string: Specifies the alignment algorithm
    fasta           // file: /path/to/genome.fasta
    gtf             // file: /path/to/genome.gtf
    bowtie2_index   // string: /path/to/bowtie2/index
    star_index      // string: /path/to/star/index

    main:

    ch_versions      = Channel.empty()
    ch_bowtie2_index = Channel.empty()
    ch_star_index    = Channel.empty()
    ch_fasta         = Channel.fromPath(fasta, checkIfExists: true).map { fasta -> [ [id: fasta.baseName], fasta ] }
    ch_gtf           = Channel.fromPath(gtf, checkIfExists: true).map { fasta -> [ [id: fasta.baseName], fasta ] }

    if (!bowtie2_index || !star_index) {

        // Append ERCC (spike-ins) to the reference genome
        spike_seq = file("$projectDir/data/spike-seq.txt")
        ERCC_CREATE(spike_seq)
        ch_versions = ch_versions.mix(ERCC_CREATE.out.versions)

        CAT_FASTA(
            ch_fasta
                .merge(ERCC_CREATE.out.fasta)
                .map{ meta, fasta, ercc -> [ [id:"${meta.id}_ercc"], [fasta, ercc] ] }
        )
        ch_versions = ch_versions.mix(CAT_FASTA.out.versions)
    }

    if (aligner.contains('bowtie2')) {
        if (bowtie2_index) {
            ch_bowtie2_index = Channel.fromPath(bowtie2_index, checkIfExists: true).map { index -> [ [id:"target_index"], index ] }
        } else {
            ch_bowtie2_index = BOWTIE2_BUILD(CAT_FASTA.out.file_out).index
            ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions)
        }
    }
    
    if (aligner.contains('star')) {
        if (star_index) {
            ch_star_index = Channel.fromPath(star_index, checkIfExists: true).map { index -> [ [id:"target_index"], index ] }
        } else {
            ch_star_index = STAR_GENOMEGENERATE( CAT_FASTA.out.file_out, ch_gtf ).index
            ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
        }
    }

    emit:
    fasta         = ch_fasta
    gtf           = ch_gtf
    bowtie2_index = ch_bowtie2_index
    star_index    = ch_star_index
    versions      = ch_versions
}