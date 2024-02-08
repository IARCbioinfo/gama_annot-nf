#! /usr/bin/env Rscript
#####################################################################################
#
# Title		: annot.r
# Author	: CahaisV@iarc.fr
# Date		: 03/07/2019
# Last Update	: 03/07/2019
#
#####################################################################################

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GetoptLong))
#suppressPackageStartupMessages(library(vcfR))
suppressPackageStartupMessages(library(Biostrings))

out=""
userAnnot=""
annovarDBpath="/data/databases/annovar/mm10db"
annovarBinPath="~/bin/annovar/"
threads=1
PASS="PASS"

GetoptLong(matrix(c("input|i=s",          "vcf input file",	
                    "annovarDBlist|l=s",  "txt file listing annovar databases for annotation",
                    "annovarDBpath|a=s",  "path to annovarDB",
                    "annovarBinPath|b=s", "path to table_annovar.pl",
                    "out|o=s",  	        "output file name",
                    "threads|t=i",        "threads",
                    "PASS|p=s",           "filter tag",
                    "keepNonStandard|k",  "keep Non Standard Chromosomes",
                    "userAnnot|u=s",      "additional user annotations"
), ncol=2, byrow=TRUE))

print(PASS)

print("Check input files")
#check annovarDBpath
if(!file.exists(annovarDBpath)){ stop( paste0("ERROR : ", annovarDBpath, " do not exit !") ) }
annovarDB=basename(annovarDBpath)
reference=gsub("db","",annovarDB)
print(paste0("Reference : ",reference ))

#check annovarDB list and create protocols and operation list
if(!file.exists(annovarDBlist)){ stop( paste0("ERROR : ", annovarDBlist, " do not exit !") ) }
dblist<-read.table(annovarDBlist,comment.char="#")
dblist$do<-apply(dblist, 1, function(row) file.exists(paste0(annovarDBpath,"/",row[1])) )
for( i in 1:nrow(dblist) ){ ifelse(dblist[i,3],print(paste0("VALID : ",dblist[i,1])),print(paste0("NOT FOUND : ", dblist[i,1])) ) }
dblist<-dblist[which(dblist$do==TRUE),]
protocols<-sub( "\\.txt", "",  sub(paste0(reference,"_"), "", dblist$V1) )
protocols<-paste(protocols,collapse = ",")
operations<-paste(dblist$V2,collapse = ",")
if(protocols==""){ stop( paste0("ERROR : there is no valid file in ", annovarDBlist, " !") ) }

#check if input file exist and read it
if(!file.exists(input)){ stop( paste0("ERROR : ", input, " do not exit !") ) }
#vcf <- read.vcfR( input, verbose = FALSE )
#colNames=c(colnames(vcf@fix), colnames(vcf@gt))
if( grepl(".gz$",input)){  
  system(command = paste("zgrep '##'", input, "> tmp") ) #save headers in tmp file
  vcf<-fread(cmd=paste("zcat", input, "| grep -v '##' "))
}else{
  system(command = paste("grep '##'", input, "> tmp") ) #save headers in tmp file
  vcf<-fread(cmd=paste("cat", input, "| grep -v '##' "))
}
colnames(vcf)[1]<-"CHROM"
colNames<-colnames(vcf)
fixNames=c("Chr","Start","End","Ref","Alt")
#keep only standard chromosomes
if (!keepNonStandard){ 
  vcf<-vcf[ !grep("_",CHROM),]
  vcf<-vcf[ !grep("M",CHROM),]

  #check genome
  nbchrom=length(unique(vcf$CHROM))
  print(paste0( "Number of chromosomes in vcf file = ", nbchrom))
  #if(reference=="hg38" & nbchrom!=24){ stop( paste0("ERROR : hg38db has been selected, but number of chromosome in vcf file is not 24 !")  ) }
  #if(reference=="hg19" & nbchrom!=24){ stop( paste0("ERROR : hg19db has been selected, but number of chromosome in vcf file is not 24 !")  ) }
  #if(reference=="mm10" & nbchrom!=21){ stop( paste0("ERROR : mm10db has been selected, but number of chromosome in vcf file is not 19 !")  ) }
  #if(reference=="mm9"  & nbchrom!=21){ stop( paste0("ERROR : mm9db has been selected,  but number of chromosome in vcf file is not 19 !")  ) }
}

#make a name for the output file
if ( out=="" ){ out=gsub("vcf.gz|vcf","tsv",input) }
print(paste0("Output file name : ", out))

###################
#filter PASS
PASS<-unlist(strsplit(PASS,","))
print(paste("PASS=",PASS))
#passvcf=gsub(".tsv","_pass.vcf.gz",out)
#write.vcf( vcf[ vcf@fix[,'FILTER'] %in% PASS ], file=passvcf )
if(!"all" %in% PASS) { vcf<-vcf[ FILTER %in% PASS] }
if(nrow(vcf)==0){ stop("ERROR : None of the variants pass filters") }
fwrite( vcf , sep="\t", file="tmp", append=TRUE)

#compress tmp vcf
err<-system("gzip -f tmp", ignore.stderr=T)
if(err>0){stop(err)}
###################
#run annovar
annovarBin<-paste0( gsub("table_annovar.pl","",annovarBinPath), "/table_annovar.pl")
params=paste(" --buildver", reference , "--thread", threads, "--vcfinput --onetranscript --remove --otherinfo --protocol", protocols, "-operation", operations, "tmp.gz", annovarDBpath , sep=" ")
print("Run table_annovar.pl")
print(paste0(annovarBin,params))
#err<-system(command, ignore.stderr=T, stdout="annovar.stdout", stderr="annovar.stderr")
err<-system2(command=annovarBin,args=params, stdout="annovar.stdout", stderr="annovar.stderr")
#if(err>0){stop(err)}
avinput<-list.files(pattern="avinput")
if(!file.exists(input)){ stop("ERROR : Annovar output not found") }
if(file.info(avinput)$size<100){ stop("ERROR : Annovar output empty") }

################
#Functions
mergeAnnovarFiles<-function(dir="./"){
  
  avinput<-list.files(path = dir, pattern="avinput")
  avinput<-fread(avinput[1])[,-c(6,7,8)]
  print(nrow(avinput))
  colnames(avinput)<-c(fixNames,colNames)
  nbfields<-ncol(avinput)
  
  variant<-list.files(path = dir, pattern="multianno.txt")
  variant<-fread(variant)
  variant<-variant[,1:(ncol(variant)-nbfields+2)]
  colnames(variant)[1:5]<-fixNames
  avinput<-merge( variant, avinput, by=fixNames, all.y=TRUE )
    
  return(avinput[ ALT!=".", ])
}


##################
#PostProcess
print("Merge annovar outputs")
avoutput<-mergeAnnovarFiles()
print(nrow(avoutput))

print("Reorder columns")
dbnames<-colnames(avoutput)[ ! colnames(avoutput) %in% c(fixNames,colNames) ]
setcolorder(avoutput, c(fixNames,dbnames,colNames) )
print(nrow(avoutput))

print("Order by chromosomic location")
avoutput$num=avoutput$Chr
avoutput[ num=="chrX", num :="23"]
avoutput[ num=="chrY", num :="24"]
avoutput[, num:=as.integer(gsub("chr","",num))]
avoutput<-avoutput[ order(c(num,Start,End,Ref,Alt)), ][!is.na(Chr)][Chr!="NA"]
avoutput[, num:=NULL]
print("Write output")
write.table(avoutput, file=out, row.names=F, sep="\t", quote=F)
print("The end !")
print(nrow(avoutput))
