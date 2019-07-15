# strelka2-nf
### Strelka v2 pipeline with Nextflow

#### Dependencies
1. Install [Strelka v2](https://github.com/Illumina/strelka).
2. Install [nextflow](http://www.nextflow.io/).

	```bash
	curl -fsSL get.nextflow.io | bash
	```
	And move it to a location in your `$PATH` (`/usr/local/bin` for example here):
	```bash
	sudo mv nextflow /usr/local/bin
	```

#### Execution

mode somatic
 `nextflow run iarcbioinfo/mutspec-annot.nf --input vcf_folder/ --annovarDBlist Dblist.txt --annovarDBpath /data/annnovar/hg38db/`

#### Options
--input 	        FOLDER	Folder containing vcf to process.
--annovarDBlist		FILE	File with two columns : protocols and operations (see example below).
--extention		TXT	input files extension
--annovarDBpath		PATH	Path to annovarDB.
--annovarBinPath	PATH	Path to table_annovar.pl.
--thread 		INT	Number of thread for table_annovar.pl.
--vaf				Add columns with VAF and coverage.


#### Help section
You can print the help manual by providing `--help` in the execution command line:
```bash
nextflow run iarcbioinfo/mutspec_annot.nf --help
```
This shows details about optional and mandatory parameters provided by the user.  

#### annovarDblist  .txt format
The annovarDBlist file is where you can define annotations. It's a tabular file with 2 columns normal and tumor.

`#This is a sample file distributed with Galaxy that is used by the
 #MutSpec-Annot tools. The mm10_listAVDB.txt has this format (white space 
 #characters are TAB characters):
 #
 #<RefGenome_DatabaseName>       <operation>
 #
 mm10_refGene.txt		g
 mm10_knownGene.txt		g
 mm10_ensGene.txt		g
 mm10_cytoBand.txt		r
 mm10_genomicSuperDups.txt	r
 mm10_snp142.txt		f`

#### Global parameters
```--annovarBinPath``` is mandatory parameters but can be defined in your nextflow config file (```~/.nextflow/config``` or ```config``` in the working directory) and so not set as inputs.

The following is an example of config part defining this:
```bash
profiles {

        standard {
                params {
                   annovarBinPath = '/data/annovar/bin/'
                }
        }
```
