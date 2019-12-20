# mutspec_annot-nf
### mutspec_annot pipeline with Nextflow

#### Dependencies
1. Install [R](https://www.r-project.org/) and libraries

	```bash
	conda install -c r r
	R
	install.packages(c("data.table","vcfR"))
	if (!requireNamespace("BiocManager", quietly = TRUE))
    		install.packages("BiocManager")
	BiocManager::install(c("biostrings","getoptlong"))
	```
	
2. Install [annovar](http://annovar.openbioinformatics.org/en/latest/user-guide/download/).
3. Install [nextflow](http://www.nextflow.io/).

	```bash
	curl -fsSL get.nextflow.io | bash
	```
	And move it to a location in your `$PATH` (`/usr/local/bin` for example here):
	```bash
	sudo mv nextflow /usr/local/bin
	```

#### Execution

 `nextflow run iarcbioinfo/mutspec-annot.nf --input vcf_folder/ --annovarDBlist Dblist.txt --annovarDBpath /data/annnovar/hg38db/`

#### Options

| OPTIONS | TYPE | Description |
|-------- | ---- | ----------- |
| --input | FOLDER | Folder containing vcf to process |
| --annovarDBlist | FILE | File with two columns : protocols and operations (see example below) |
| --extention | TXT | input files extension |
| --annovarDBpath | PATH | Path to annovarDB |
| --annovarBinPath | PATH | Path to table_annovar.pl |
| --thread | INT | Number of thread for table_annovar.pl |
| --vaf |  | Add columns with VAF and coverage |
| --caller | TXT | when using --vaf, indicate the caller (strelka2, mutect2 or haplotypecaller) |

#### Help section
You can print the help manual by providing `--help` in the execution command line:
```bash
nextflow run iarcbioinfo/mutspec_annot.nf --help
```
This shows details about optional and mandatory parameters provided by the user.  

#### annovarDblist  .txt format
The annovarDBlist file is where you can define annotations. See example [mm10.list](https://github.com/IARCbioinfo/mutspec_annot/blob/master/demo/mm10.list) in demo folder.

### annovarDBpath 
The annovarDBpath is where your annovar database is located. This folder is created using the annotate_variation.pl script from annovar. Plus you need to add int the same folder the fasta of your genome for the context annotation.

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
