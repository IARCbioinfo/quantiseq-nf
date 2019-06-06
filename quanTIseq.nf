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
params.output_folder= "."
params.mem  = 2
params.cpu  = 1
params.suffix1   = "_1"
params.suffix2   = "_2"
params.fastq_ext = "fastq.gz"
params.image = null

params.help = null

log.info ""
log.info "--------------------------------------------------------"
log.info "  quantiseq-nf <VERSION>: <SHORT DESCRIPTION>         "
log.info "--------------------------------------------------------"
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
    log.info '    --input_folder   FOLDER                  Folder containing fastq files.'
    log.info '    --suffix1        STRING                  Suffix for fastq file with 1st element of pair.'
    log.info '    --suffix2        STRING                  Suffix for fastq file with 2nd element of pair.'
    log.info ""
    log.info "Optional arguments:"
    log.info '    --output_folder     STRING                Output folder (default: .).'
    log.info '    --cpu          INTEGER                 Number of cpu used (default: 1).'
    log.info '    --mem          INTEGER                 Size of memory (in GB) (default: 2).' 
    exit 0
} else {
/* Software information */
   log.info "input_folder = ${params.input_folder}"
   log.info "cpu          = ${params.cpu}"
   log.info "mem          = ${params.mem}"
   log.info "output_folder= ${params.output_folder}"
   log.info "help:                               ${params.help}"
}

if (file(params.input_folder).listFiles().findAll { it.name ==~ /.*${params.fastq_ext}/ }.size() > 0){
    println "fastq files found, proceed with quantification"
}else{
	println "ERROR: input folder contains no fastq files"; System.exit(0)
}

keys1 = file(params.input_folder).listFiles().findAll { it.name ==~ /.*${params.suffix1}.${params.fastq_ext}/ }.collect { it.getName() }
                                                                                                               .collect { it.replace("${params.suffix1}.${params.fastq_ext}",'') }
keys2 = file(params.input_folder).listFiles().findAll { it.name ==~ /.*${params.suffix2}.${params.fastq_ext}/ }.collect { it.getName() }
                                                                                                               .collect { it.replace("${params.suffix2}.${params.fastq_ext}",'') }
if ( !(keys1.containsAll(keys2)) || !(keys2.containsAll(keys1)) ) {println "\n ERROR : There is not at least one fastq without its mate, please check your fastq files."; System.exit(0)}
println keys1

// Gather files ending with _1 suffix
reads1 = Channel
    .fromPath( params.input_folder+'/*'+params.suffix1+'.'+params.fastq_ext )
    .map {  path -> [ path.name.replace("${params.suffix1}.${params.fastq_ext}",""), path ] }

// Gather files ending with _2 suffix
reads2 = Channel
    .fromPath( params.input_folder+'/*'+params.suffix2+'.'+params.fastq_ext )
    .map {  path -> [ path.name.replace("${params.suffix2}.${params.fastq_ext}",""), path ] }

// Match the pairs on two channels having the same 'key' (name) and emit a new pair containing the expected files
readPairs = reads1
    .phase(reads2)
    .map { pair1, pair2 -> [ pair1[1], pair2[1] ] }

    println reads1


if(params.image!=null){
	image=file(params.image)
}else{
	process pullsingularity {
	        cpus 1
	        memory '1G'
	        tag { file_tag }

	        output:
	        file("quantiseq2.img") into image

	        publishDir "${params.output_folder}", mode: 'copy'

	        shell:
	        '''
		singularity pull IARCbioinfo/quantiseq-nf:v4
		mv  IARCbioinfo-quantiseq-nf-master-v4.simg quantiseq2.img
	        '''
	}
}

// launches quanTIseq
process quanTIseq {
	cpus params.cpu
	memory params.mem+'G'
	tag { file_tag }

	input:
	file pairs from readPairs
	file image
	
	output:
	file("*.txt") into outputs
	file("*.tsv") into outputs2

	publishDir "${params.output_folder}", mode: 'copy'
    
	shell:
	file_tag = pairs[0].name.replace("${params.suffix1}.${params.fastq_ext}","")
    	'''
	echo "!{file_tag}\t!{pairs[0]}\t!{pairs[1]}" > input.txt
	!{baseDir}/bin/quanTIseq_pipeline.sh --prefix=quanTIseq_!{file_tag} --threads=!{params.cpu} --inputfile=input.txt --outputdir=!{params.input}
    	'''
}
