params.help = null
params.input = null
params.extension = "vcf"
params.annovarDBlist  = null
params.annovarDBpath  = "/data/databases/annovar/hg38db/"
params.annovarBinPath = "/data/mca/mca_share/work/annovar/"
params.output = "mutspec_annotation"
params.thread = 1
params.VAF = true
params.caller = "strelka2"
#params.pass = "'PASS,clustered_events,clustered_events;homologous_mapping_event,tiers1,tiers2,tiers3'"
params.pass = "'PASS'"

if (params.help) {
    log.info ''
    log.info '--------------------------------------------------'
    log.info '                   TABLE ANNOVAR                  '
    log.info '--------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info 'nextflow run mutspec_annot.nf --input vcfFolder/ --annovarDBpath /data/annovar/mm10db/'
    log.info ''
    log.info 'Mandatory arguments:'
    log.info ''
    log.info '    --input            FOLDER            Folder containing vcf to process.'
    log.info '    --annovarDBlist    FILE              File with two columns : protocols and operations.'
    log.info ''
    log.info 'Optional arguments:'
    log.info ''
    log.info '    --extention        STRING            input files extension'
    log.info '    --annovarDBpath    PATH              Path to annovarDB.'
    log.info '    --annovarBinPath   PATH              Path to table_annovar.pl.'
    log.info '    --output           FOLDER            Output Folder name.'
    log.info '    --thread           INTEGER           Number of thread for table_annovar.pl.'
    log.info '    --caller           PATH              Software used for calling (strelka2, mutect2 or haplotypecaller)'
    log.info '    --pass             STRING            filter tags, comma separated list'
    log.info ''   
    log.info 'Flags'
    log.info ''
    log.info '   --vaf                                 Add columns with VAF and coverage.'
    log.info '   --help                                Display this message.'
    exit 1
}



//vcf=Channel.fromFilePairs( params.input + '/*{snvs,indels}*' + params.extension )
//System.exit(0)
tsize=2
Channel.fromPath( params.input + '/*indels*' + params.extension ).ifEmpty { tsize=1 }
allvcf = Channel.fromPath( params.input + '/*' + params.extension ).ifEmpty { error "empty table folder, please verify your input." }


process mutspec_annot {

  publishDir params.output, mode: 'copy'

  cpus params.thread
  tag { sample_tag }

  input:
  file vcf from allvcf

  output:
  set val(sample_tag), file("*tab") into annotated

  shell:
  sample_tag = vcf.baseName.replaceFirst(/(snvs|indels).*/,"")
  '''
  echo !{tsize}
  echo !{sample_tag}
  echo mutspec_annot.r -i !{vcf} -l !{params.annovarDBlist} -a !{params.annovarDBpath} -b !{params.annovarBinPath} -t !{params.thread} -p "!{params.pass}"
  mutspec_annot.r -i !{vcf} -l !{params.annovarDBlist} -a !{params.annovarDBpath} -b !{params.annovarBinPath} -t !{params.thread} -p "!{params.pass}"
  '''

}

if (params.VAF){

  process mutspec_VAF {

    publishDir params.output, mode: 'move'

    input:
    set val(sample_tag), file(tab) from annotated.groupTuple(size: tsize )

    output:
    file "*tab" into table

    shell:
    '''
    getAllelicFraction.r -c !{params.caller}
    '''
  }

}






