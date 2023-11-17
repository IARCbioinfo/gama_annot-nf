#! /usr/bin/env nextflow

// Copyright (C) 2017 IARC/WHO

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

params.help = null
params.input_folder = null
params.annovarDBlist  = null
params.annovarDBpath  = "/data/databases/annovar/hg38db/"
params.annovarBinPath = "~/bin/annovar/"
params.output_folder = "gama_annot"
params.cpu = 8
params.mem = 64
params.pass = "'PASS'"
//params.pass = "'PASS,clustered_events,clustered_events;homologous_mapping_event,tiers1,tiers2,tiers3'"


log.info ""
log.info "------------------------------------------------------------------------"
log.info "  Gama-annot 1.2 : Annotation of vcf files with annovar (nextflow DSL2) "
log.info "------------------------------------------------------------------------"
log.info "Copyright (C) IARC/WHO"
log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
log.info "This is free software, and you are welcome to redistribute it"
log.info "under certain conditions; see LICENSE for details."
log.info "--------------------------------------------------------"
log.info ""

if (params.help) {
    log.info ''
    log.info 'Usage: '
    log.info 'nextflow run gama_annot.nf --input vcfFolder/ --annovarDBlist dblist --annovarDBpath /data/annovar/hg38db/'
    log.info ''
    log.info 'Mandatory arguments:'
    log.info ''
    log.info '    --input            FOLDER            Folder containing vcf to process.'
    log.info '    --annovarDBlist    FILE              File with two columns : protocols and operations (see anovar documentation).'
    log.info ''
    log.info 'Optional arguments:'
    log.info ''
    log.info '    --annovarDBpath    PATH              Path to annovarDB.'
    log.info '    --annovarBinPath   PATH              Path to table_annovar.pl.'
    log.info '    --output           FOLDER            Output Folder name.'
    log.info '    --cpu              INTEGER           Number of cpu used by table_annovar.pl default (8)'
    log.info '    --mem              INTEGER           Size of memory used by gama_annot in GB default (64)'
    log.info '    --pass             STRING            filter tags, comma separated list'
    log.info ''   
    log.info 'Flags'
    log.info ''
    log.info '   --help                                Display this message.'
    exit 1
}

/***************************************************************************************/
/************************ handle global parameters *************************************/
/***************************************************************************************/


/***************************************************************************************/
/************************  Process   ***************************************************/
/***************************************************************************************/

process gama_annot {

  memory = params.mem+'.GB'
  cpus params.cpu

  input:
    path vcf

  output:
    tuple val(sample_tag), path("*tsv"), emit: annotated

  shell:
    sample_tag = vcf.baseName.replaceFirst(/(.snvs|.indels).*/,"")
    """
    annot.r -i ${vcf} -t ${params.cpu} -p "${params.pass}" \
            -l ${params.annovarDBlist} -a ${params.annovarDBpath} -b ${params.annovarBinPath}
    """

  stub:
    sample_tag = vcf.baseName.replaceFirst(/(.snvs|.indels).*/,"")
    """
    touch ${sample_tag}.tsv
    """

}


process gama_context {
    
    publishDir params.output_folder, mode: 'copy'

    memory = params.mem+'.GB'
    cpus params.cpu

    input:
      tuple val(sample_tag), path(tab)

    output:
      tuple val(sample_tag), file("*1.tsv"), emit: context

    shell:
      """
      getContext.r -a ${params.annovarDBpath}
      """

    stub:
      """
      touch ${sample_tag}.1.tsv
      """
}



/****************************************************************************************/
/************************  Workflow   ***************************************************/
/****************************************************************************************/

workflow {

  allvcf = Channel.fromPath( "${params.input_folder}/*.{vcf,vcf.gz}" ).ifEmpty { 
    error "empty table folder, please verify your input." 
  }

  gama_annot(allvcf) | gama_context

}