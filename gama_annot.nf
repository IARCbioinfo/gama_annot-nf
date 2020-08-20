params.help = null
params.input = null
params.extension = "vcf"
params.annovarDBlist  = null
params.annovarDBpath  = "/data/databases/annovar/hg38db/"
params.annovarBinPath = "~/bin/annovar/"
params.output = "gama_annot"
params.cpu = 8
params.mem = 64
params.context = false
params.caller = "none"
params.pass = "'PASS'"
//params.pass = "'PASS,clustered_events,clustered_events;homologous_mapping_event,tiers1,tiers2,tiers3'"

if (params.help) {
    log.info ''
    log.info '--------------------------------------------------'
    log.info '                   TABLE ANNOVAR                  '
    log.info '--------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info 'nextflow run gama_annot.nf --input vcfFolder/ --annovarDBpath /data/annovar/mm10db/'
    log.info ''
    log.info 'Mandatory arguments:'
    log.info ''
    log.info '    --input            FOLDER            Folder containing vcf to process.'
    log.info '    --annovarDBlist    FILE              File with two columns : protocols and operations (see anovar documentation).'
    log.info ''
    log.info 'Optional arguments:'
    log.info ''
    log.info '    --extention        STRING            input files extension'
    log.info '    --annovarDBpath    PATH              Path to annovarDB.'
    log.info '    --annovarBinPath   PATH              Path to table_annovar.pl.'
    log.info '    --output           FOLDER            Output Folder name.'
    log.info '    --cpu              INTEGER           Number of cpu used by table_annovar.pl default (8)'
    log.info '    --mem              INTEGER           Size of memory used by gama_annot in GB default (64)'
    log.info '    --caller           STRING            Add columns with VAF and coverage for corresponding software (strelka2, mutect2 or haplotypecaller)'
    log.info '    --pass             STRING            filter tags, comma separated list'
    log.info ''   
    log.info 'Flags'
    log.info ''
    log.info '   --context                             Add context, strand and user Annotations'
    log.info '   --help                                Display this message.'
    exit 1
}


tsize=2
Channel.fromPath( params.input + '/*indels*' + params.extension ).ifEmpty { tsize=1 }
allvcf = Channel.fromPath( params.input + '/*' + params.extension ).ifEmpty { error "empty table folder, please verify your input." }


process gama_annot {

  publishDir params.output, mode: 'copy'

  memory = params.mem+'.GB'
  cpus params.cpu
  tag { sample_tag }

  input:
  file vcf from allvcf

  output:
  set val(sample_tag), file("*tsv") into annotated

  shell:
  sample_tag = vcf.baseName.replaceFirst(/(snvs|indels).*/,"")
  '''
  echo !{tsize}
  echo !{sample_tag}
  echo annot.r -i !{vcf} -l !{params.annovarDBlist} -a !{params.annovarDBpath} -b !{params.annovarBinPath} -t !{params.cpu} -p "!{params.pass}"
  annot.r -i !{vcf} -l !{params.annovarDBlist} -a !{params.annovarDBpath} -b !{params.annovarBinPath} -t !{params.cpu} -p "!{params.pass}"
  '''

}

if (params.context){

process gama_context {
    
    publishDir params.output, mode: 'copy'

    memory = params.mem+'.GB'
    cpus params.cpu

    input:
    set val(sample_tag), file(tab) from annotated

    output:
    set val(sample_tag), file("*1.tsv") into res_context

    shell:
    '''
    getContext.r -a !{params.annovarDBpath}
    '''
}

}else{
   res_context = annotated
}

if (params.caller!="none"){

  process gama_VAF {

    publishDir params.output, mode: 'move'

    memory = params.mem+'.GB'
    cpus params.cpu
    
    input:
    set val(sample_tag), file(tab) from res_context.groupTuple(size: tsize )

    output:
    file "*2.tsv" into table

    shell:
    '''
    getAllelicFraction.r -c !{params.caller}
    '''
  }

}






