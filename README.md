# gama_annot-nf
### gama_annot pipeline with Nextflow

Annotate vcf files using annovar + additional scripts. Works for Mutect2, Strelka2 and HaplotypeCaller outputs

#### Dependencies
1. Install [R](https://www.r-project.org/) and libraries

	```bash
	conda install -c r r
	R
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

4. If the vcf files to annotate are from Strelka2 but were not generated using IARCbioinfo Strelka2 pipeline, first run  [fixStrelkaOutput.sh](https://github.com/IARCbioinfo/strelka2-nf/blob/master/bin/fixStrelkaOutput.sh). This will add GT field (mandatory for Annovar)

#### Execution

 `nextflow run iarcbioinfo/gama_annot-nf --input vcf_folder/ --annovarDBlist Dblist.txt --annovarDBpath /data/annnovar/hg38db/ --pass PASS --context --caller strelka2 --extention vcf.gz

  nextflow run iarcbioinfo/gama_annot-nf --input vcf_folder/ --annovarDBlist Dblist.txt --annovarDBpath /data/annnovar/hg38db/ --pass PASS,clustered_events --context --caller mutect2 --extention vcf

#### Options

| OPTIONS | TYPE | Description |
|-------- | ---- | ----------- |
| --input | FOLDER | Folder containing vcf to process |
| --annovarDBlist | FILE | File with two columns : protocols and operations (see example below) |
| --extension | TXT | input files extension |
| --annovarDBpath | PATH | Path to annovarDB |
| --annovarBinPath | PATH | Path to table_annovar.pl |
| --thread | INT | Number of thread for table_annovar.pl |
| --caller | TXT | when using --vaf, indicate the caller (strelka2, mutect2 or haplotypecaller) |
| --pass | TXT | Value on which the variants should be filtered prior to annotation (default : PASS) | 

#### Help section
You can print the help manual by providing `--help` in the execution command line:
```bash
nextflow run iarcbioinfo/gama_annot-nf --help
```
This shows details about optional and mandatory parameters provided by the user.  

#### annovarDblist  .txt format
The annovarDBlist file is where you can define annotations. See example [mm10.list](https://github.com/IARCbioinfo/gama_annot-nf/blob/master/demo/mm10.list) in demo folder.

### annovarDBpath 
The annovarDBpath is where your annovar database is located. This folder is created using the annotate_variation.pl script from annovar. The name of the folder should be genomedb (for example hg38db or mm10db).
You also need to add in the same folder the reference of your genome for the context annotation (specific ".fa" file compatible with the specific context annotation)

#### Global parameters
```--annovarBinPath``` is mandatory parameters but can be defined in your nextflow config file (```~/.nextflow/config``` or ```config``` in the working directory) and so not set as inputs.

The following is an example of config part defining this:
```bash
profiles {

        standard {
                params {
                   annovarBinPath = '/data/annovar/bin/'
                   annovarDBpath = '/data/annovarDB/hg38db/'
                }
        }
```
