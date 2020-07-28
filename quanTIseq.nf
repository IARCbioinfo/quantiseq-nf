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
log.info "  quantiseq-nf v1.1: quantification of immune infiltration with software quanTIseq "
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

if(params.input_file){
	readPairs_premerge = Channel.fromPath("${params.input_file}")
     	   .splitCsv( header: true, sep: '\t', strip: true )
	       .map { row -> [ row.SM , file(row.pair1), file(row.pair2) ] }
           .groupTuple(by: 0)
}else{
    if (file(params.input_folder).listFiles().findAll { it.name ==~ /.*${params.fastq_ext}/ }.size() > 0){
        println "fastq files found, proceed with quantification"
    }else{
	    println "ERROR: input folder contains no fastq files"; System.exit(0)
    }

    readPairs = Channel.fromFilePairs(params.input_folder +"/*{${params.suffix1},${params.suffix2}}" +'.'+ params.fastq_ext)
			      .map { row -> [ row[0] , row[1][0], row[1][1] ] }
}

if(params.image!=null){
	image=file(params.image)
}else{
	process pullsingularity {
	        cpus 1
	        memory '1G'

	        output:
	        file("quantiseq2.img") into image

	        publishDir "${params.output_folder}", mode: 'copy'

	        shell:
	        '''
		    singularity pull quantiseq2.img shub://IARCbioinfo/quantiseq-nf:v1.1
	        '''
	}
}

if(params.input_file){
    readPairsNot2merge = Channel.create()
    readPairs2merge    = Channel.create()
    readPairs_premerge.choice( readPairsNot2merge, readPairs2merge ) { a -> a[1].size()==1 ? 0 : 1 }

    process merge {
	        cpus 2
	        memory params.mem+'G'
	        tag { SM }

            input:
	        set val(SM), file(pair1), file(pair2) from readPairs2merge

	        output:
	        set val(SM), file("${SM}${params.suffix1}.${params.fastq_ext}"), file("${SM}${params.suffix2}.${params.fastq_ext}") into readPairsMerged

	        shell:
	        '''
		    cat !{pair1} > !{SM}!{params.suffix1}.!{params.fastq_ext}
            cat !{pair2} > !{SM}!{params.suffix2}.!{params.fastq_ext}
	        '''
	}
    readPairs = readPairsNot2merge.concat( readPairsMerged )
}

// launches quanTIseq
process quanTIseq {
	cpus params.cpu
	memory params.mem+'G'
	tag { SM }

	input:
	set val(SM), file(pair1), file(pair2) from readPairs
	file image
	
	output:
	file("quantiseqResults*/*txt") into res_quantiseq

	publishDir "${params.output_folder}/intermediate_results/", mode: 'copy'
    
	shell:
	mode = params.nontumor  ? "" : "--tumor=TRUE"
    '''
	echo "!{SM}\t!{pair1}\t!{pair2}" > input.txt
	!{baseDir}/bin/quanTIseq_pipeline.sh --threads=!{params.cpu} --inputfile=input.txt --outputdir=. !{mode}
    mv quantiseqResults_* quantiseqResults_!{SM} 
    cd quantiseqResults_!{SM}/ && mv quanTIseq_cell_fractions.txt quanTIseq_cell_fractions_!{SM}.txt && mv quanTIseq_gene_tpm.txt quanTIseq_gene_tpm_!{SM}.txt
    '''
}

process merge_quanTIseq_res {
	cpus 2
	memory '300M'

	input:
	file res from res_quantiseq.collect()
	
	output:
	file("quanTIseq_*matrix.txt") into output_mat

	publishDir "${params.output_folder}", mode: 'copy'
    
	shell:
    '''
    awk 'FNR==1 && NR!=1 { while (/^Sample/) getline; } 1 {print}' quanTIseq_cell_fractions*.txt > quanTIseq_cell_fractions_matrix.txt
    for f in `ls quanTIseq_gene_tpm_*.txt`; do cut -f2 $f > $f.cut && cut -f1 $f > rownames_gene_tpm.txt; done
    paste rownames_gene_tpm.txt quanTIseq_gene_tpm_*.txt.cut > quanTIseq_gene_tpm_matrix.txt
	'''
}