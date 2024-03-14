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
params.tags = "'.bam|.vcf|.tsv|.somatic|.germline|.calls|.alt|.gz|.PASS|.filtered'"
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
    log.info '    --input_folder     FOLDER            Folder containing vcf to process.'
    log.info '    --annovarDBlist    FILE              File with two columns : protocols and operations (see anovar documentation).'
    log.info ''
    log.info 'Optional arguments:'
    log.info ''
    log.info '    --annovarDBpath    PATH              Path to annovarDB.'
    log.info '    --annovarBinPath   PATH              Path to table_annovar.pl.'
    log.info '    --output_folder    FOLDER            Output Folder name.'
    log.info '    --cpu              INTEGER           Number of cpu used by table_annovar.pl default (8)'
    log.info '    --mem              INTEGER           Size of memory used by gama_annot in GB default (64)'
    log.info '    --pass             STRING            filter flags, comma separated list'
    log.info '    --tags             STRING            strings paterns to remove in input file name'
    log.info ''   
    log.info 'Flags'
    log.info ''
    log.info '   --help                                Display this message.'
    exit 1
}

/***************************************************************************************/
/************************ handle global parameters *************************************/
/***************************************************************************************/

log.info ""
log.info "input_folder          = ${params.input_folder}"
log.info "annovarDBlist         = ${params.annovarDBlist}"
log.info "annovarDBpath         = ${params.annovarDBpath}"
log.info "annovarBinPath        = ${params.annovarBinPath}"
log.info "output_folder         = ${params.output_folder}"
log.info "cpu                   = ${params.cpu}"
log.info "mem                   = ${params.mem}"
log.info "pass                  = ${params.pass}"
log.info "tags                  = ${params.tags}"
log.info ""

/***************************************************************************************/
/************************  Process   ***************************************************/
/***************************************************************************************/

process annovar_annot {

  memory = params.mem+'.GB'
  cpus params.cpu

  publishDir params.output_folder, mode: 'copy', pattern: '{*multianno*}'

  input:
    path vcf
    path annovar
    path annovarDB

  output:
    tuple val(sample_tag), path("*avinput"), path("*multianno.txt"), path("*multianno.vcf"), emit: annotated

  shell:
    sample_tag = vcf.baseName.replaceFirst(/(.snvs|.indels).*/,"")
    """
    annovar_annot.r -i ${vcf} -t ${params.cpu} -p "${params.pass}" -g "${params.tags}" \
            -l ${params.annovarDBlist} -a ${annovarDB} -b ${annovar}
    """

  stub:
    sample_tag = vcf.baseName.replaceFirst(/(.snvs|.indels).*/,"")
    """
    touch ${sample_tag}.tsv
    """

}


process gama_annot {
    
    publishDir params.output_folder, mode: 'copy'

    memory = params.mem+'.GB'
    cpus params.cpu

    input:
      tuple val(sample_tag), path(avinput), path(tab), path(vcf)
      path annovarDB

    output:
      tuple val(sample_tag), file("*1.tsv"), emit: context

    shell:
      """
      gama_annot.r -a ${annovarDB}
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

  def annovar = file( params.annovarBinPath )
  def annovarDB = file( params.annovarDBpath )

  allvcf = Channel.fromPath( "${params.input_folder}/*.{vcf,vcf.gz}" ).ifEmpty { 
    error "empty table folder, please verify your input." 
  }

  annovar_annot(allvcf,annovar,annovarDB)
  gama_annot(annovar_annot.out.annotated,annovarDB)

}
