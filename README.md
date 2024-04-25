# gama_annot-nf

## gama_annot pipeline with Nextflow

Annotate vcf files using annovar + additional scripts. Works for Mutect2, Strelka2 and HaplotypeCaller outputs

### Dependencies

1. Nextflow : for common installation procedures see the [IARC-nf](https://github.com/IARCbioinfo/IARC-nf) repository.

2. Install [annovar](http://annovar.openbioinformatics.org/en/latest/user-guide/download/).

3. Environment

	**A conda receipe, and docker and singularity containers are available with all the tools needed to run the pipeline (see "Usage")**

	```bash
	conda env create -f environment.yml
	```

### Use Annovar with strelka vcf

If the vcf files to annotate are from Strelka2 but were not generated using IARCbioinfo Strelka2 pipeline, first run  [fixStrelkaOutput.sh](https://github.com/IARCbioinfo/strelka2-nf/blob/master/bin/fixStrelkaOutput.sh). This will add GT field (mandatory for Annovar)

### Execution

 `nextflow run iarcbioinfo/gama_annot-nf -r master -latest -profile singularity --annovarDBlist Dblist.txt --annovarDBpath /data/annnovar/hg38db/ --annovarBinPath ~/bin/annovar/ --input_folder vcf_folder/`

### Options

| OPTIONS | TYPE   | Description |
|-------- | ------ | ----------- |
| --input_folder | FOLDER | Folder containing vcf to process |
| --annovarDBlist  | FILE | File with two columns : protocols and operations (see example below) |
| --annovarDBpath  | PATH | Path to annovarDB |
| --annovarBinPath | PATH | Path to table_annovar.pl |
| --pass | STRING | filter flags, comma separated list |
| --tags | STRING | tags to remove in input file names |
| --cpu  | INT | Number of used by table_annovar.pl default (8) |
| --mem  | INT | Size of memory used by gama_annot in GB default (64) |



### Help section
You can print the help manual by providing `--help` in the execution command line:

```bash
nextflow run iarcbioinfo/gama_annot-nf --help
```
This shows details about optional and mandatory parameters provided by the user.  


### annovarBinPath

This is the location of annovar perl scripts on your system.

### annovarDBpath 

The annovarDBpath is where your annovar database is located. This folder is created using the annotate_variation.pl script from annovar. The name of the folder should be genomedb (for example hg38db or mm10db).
You also need to add in the same folder the reference of your genome for the context annotation (specific ".fa" file compatible with the specific context annotation)

### annovarDblist

The annovarDBlist file is where you can choose annotations databases. See example [hg38_listAVDB.txt](https://github.com/IARCbioinfo/gama_annot-nf/blob/master/demo/hg38_listAVDB.txt) in demo folder. Each line is database accessible in annovarDBpath.


### profile

The following is an example of config for gama_annot:

```bash
profiles {

        hg38 {
                params {
                   annovarBinPath = '/data/annovar/bin/'
                   annovarDBpath = '/data/annovarDB/hg38db/'
				   annovarDBlist = '/data/annovar/hg38_listAVDB.txt'
                }
        }		
```
