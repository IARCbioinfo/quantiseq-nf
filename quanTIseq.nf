#! /usr/bin/env nextflow

//vim: syntax=groovy -*- mode: groovy;-*-

// Copyright (C) 2018 IARC/WHO

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

nextflow.enable.dsl = 2

// --------------------------------------------------
// PARAMETERS
// --------------------------------------------------

params.input_folder = "."
params.input_file   = null
params.output_folder= "."
params.mem  = 2
params.cpu  = 1
params.suffix1   = "_1"
params.suffix2   = "_2"
params.fastq_ext = "fq.gz"
params.image = null
params.nontumor = null
params.help = null

log.info ""
log.info "-----------------------------------------------------------------------------------"
log.info "  quantiseq-nf v1.1: quantification of immune infiltration with quanTIseq"
log.info "-----------------------------------------------------------------------------------"
log.info "Copyright (C) IARC/WHO"
log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
log.info "This is free software, and you are welcome to redistribute it"
log.info "under certain conditions; see LICENSE for details."
log.info "--------------------------------------------------------"
log.info ""

if (params.help) {
    log.info "--------------------------------------------------------"
    log.info "  USAGE                                                 "
    log.info "--------------------------------------------------------"
    log.info ""
    log.info "nextflow run iarcbioinfo/rnaseq-transcript-nf [-with-docker] [OPTIONS]"
    log.info ""
    log.info "Mandatory arguments:"
    log.info '    --input_folder   FOLDER              Folder containing fastq files.'
    log.info ""
    log.info "Optional arguments:"
    log.info '    --input_file      STRING             Input file (tab-separated values) with 3 columns:'
    log.info '                                         SM (sample name), pair1 (first fastq pair file),'
    log.info '                                         and pair2 (second fastq pair file).'
    log.info '    --output_folder   STRING             Output folder (default: .).'
    log.info '    --suffix1         STRING             Suffix for fastq file with 1st element of pair.'
    log.info '    --suffix2         STRING             Suffix for fastq file with 2nd element of pair.'
    log.info '    --fastq_ext       STRING             Extension of fastq files (default : fq.gz)'
    log.info '    --cpu             INTEGER            Number of cpu used (default: 1).'
    log.info '    --mem             INTEGER            Size of memory (in GB) (default: 2).' 
    log.info '    --image           STRING             Path to quantiseq singularity image (default: null).' 
    log.info ''
    log.info 'Flags:'
    log.info '    --nontumor                           Use nontumor quantiseq mode'
    exit 0
} else {
/* Software information */
   log.info "input_folder = ${params.input_folder}"
   log.info "input_file     = ${params.input_file}"
   log.info "cpu          = ${params.cpu}"
   log.info "mem          = ${params.mem}"
   log.info "suffix1      = ${params.suffix1}"
   log.info "suffix2      = ${params.suffix2}"
   log.info "output_folder= ${params.output_folder}"
   log.info "fastq_ext    = ${params.fastq_ext}"
   log.info "image        = ${params.image}"
   log.info "nontumor     = ${params.nontumor}"
   log.info "help:        ${params.help}"
}

// --------------------------------------------------
// PROCESSES
// --------------------------------------------------

process PULLSINGULARITY {
    cpus 1
    memory '1G'

    output:
    path "quantiseq2.img"

    publishDir "${params.output_folder}", mode: 'copy'

    script:
    '''
    singularity pull quantiseq2.img shub://IARCbioinfo/quantiseq-nf:v1.1
    '''
}

process MERGE_FASTQ {
    tag { SM }
    cpus 2
    memory "${params.mem}G"

    input:
    tuple val(SM), path(pair1), path(pair2)

    output:
    tuple val(SM),
          path("${SM}${params.suffix1}.${params.fastq_ext}"),
          path("${SM}${params.suffix2}.${params.fastq_ext}")

    script:
    '''
    cat ${pair1} > ${SM}${params.suffix1}.${params.fastq_ext}
    cat ${pair2} > ${SM}${params.suffix2}.${params.fastq_ext}
    '''
}

process QUANTISEQ {
    tag { SM }
    cpus params.cpu
    memory "${params.mem}G"

    input:
    tuple val(SM), path(pair1), path(pair2)
    path image

    output:
    path "quantiseqResults*/*txt"

    publishDir "${params.output_folder}/intermediate_results", mode: 'copy'

    script:
    def mode = params.nontumor ? "" : "--tumor=TRUE"
    """
    echo "${SM}\t${pair1}\t${pair2}" > input.txt
    ${baseDir}/bin/quanTIseq_pipeline.sh \
        --threads=${params.cpu} \
        --inputfile=input.txt \
        --outputdir=. ${mode}

    mv quantiseqResults_* quantiseqResults_${SM}
    cd quantiseqResults_${SM}
    mv quanTIseq_cell_fractions.txt quanTIseq_cell_fractions_${SM}.txt
    mv quanTIseq_gene_tpm.txt quanTIseq_gene_tpm_${SM}.txt
    """
}

process MERGE_QUANTISEQ_RESULTS {
    cpus 2
    memory '300M'

    input:
    path res_files

    output:
    path "quanTIseq_*matrix.txt"

    publishDir "${params.output_folder}", mode: 'copy'

    script:
    '''
    awk 'FNR==1 && NR!=1 { while (/^Sample/) getline; } 1 {print}' \
        quanTIseq_cell_fractions*.txt > quanTIseq_cell_fractions_matrix.txt

    for f in quanTIseq_gene_tpm_*.txt; do
        cut -f2 $f > $f.cut
        cut -f1 $f > rownames_gene_tpm.txt
    done

    paste rownames_gene_tpm.txt quanTIseq_gene_tpm_*.txt.cut \
        > quanTIseq_gene_tpm_matrix.txt
    '''
}

// --------------------------------------------------
// WORKFLOW
// --------------------------------------------------

workflow {

    // 1) Input ///////////////////////////////////////

    def readPairs_premerge
    def readPairs

    if (params.input_file) {
        readPairs_premerge =
            Channel.fromPath(params.input_file)
                   .splitCsv(header: true, sep: '\t', strip: true)
                   .map { row -> tuple(row.SM, file(row.pair1), file(row.pair2)) }
                   .groupTuple(by: 0)
    } else {
        readPairs =
            Channel.fromFilePairs(
                "${params.input_folder}/*{${params.suffix1},${params.suffix2}}.${params.fastq_ext}"
            ).map { row ->
                tuple(row[0], row[1][0], row[1][1])
            }
    }

    // 2) Singularity image ////////////////////////////

    def image_ch

    if (params.image) {
        image_ch = Channel.value(file(params.image))
    } else {
        image_ch = PULLSINGULARITY.out
    }

    // 3) Merge FASTQs if needed ///////////////////////

    if (params.input_file) {

        def readPairsNot2merge
        def readPairs2merge

        readPairs_premerge
            .choice(readPairsNot2merge, readPairs2merge) { it[1].size() == 1 ? 0 : 1 }

        def merged = MERGE_FASTQ(readPairs2merge)

        readPairs = readPairsNot2merge.concat(merged)
    }

    // 4) Run quanTIseq ////////////////////////////////
	
    def quant_results = QUANTISEQ(readPairs, image_ch)

    // 5) Merge outputs ////////////////////////////////

    MERGE_QUANTISEQ_RESULTS(quant_results.collect())
}
