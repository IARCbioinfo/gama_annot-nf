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
suppressPackageStartupMessages(library(vcfR))
suppressPackageStartupMessages(library(Biostrings))

out=""
annovarDBpath="/data/databases/annovar/mm10db"
annovarBinPath="~/bin/"
threads=1
PASS="PASS"

GetoptLong(matrix(c("input|i=s",          "vcf input file",	
                    "annovarDBlist|l=s",  "txt file listing annovar databases for annotation",
                    "annovarDBpath|a=s",  "path to annovarDB",
                    "annovarBinPath|b=s", "path to table_annovar.pl",
                    "out|o=s",  	  "output file name",
                    "threads|t=s",        "threads",
                    "PASS|p=s",           "filter tag"
), ncol=2, byrow=TRUE))

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
vcf <- read.vcfR( input, verbose = FALSE )
colNames=c(colnames(vcf@fix), colnames(vcf@gt))
fixNames=c("Chr","Start","End","Ref","Alt")

#check genome
nbchrom=length(unique(vcfR::getCHROM(vcf)))
print(paste0( "Number of chromosomes in vcf file = ", nbchrom))
if(reference=="hg38" & nbchrom!=24){ print( paste0("WARNING : hg38db has been selected, but number of chromosome is not 24 !")  ) }
if(reference=="hg19" & nbchrom!=24){ print( paste0("WARNING : hg19db has been selected, but number of chromosome is not 24 !")  ) }
if(reference=="mm10" & nbchrom!=21){ print( paste0("WARNING : mm10db has been selected, but number of chromosome is not 19 !")  ) }
if(reference=="mm9"  & nbchrom!=21){ print( paste0("WARNING : mm9db has been selected, but number of chromosome is not 19 !")  ) }

#make a name for the output file
if ( out=="" ){ out=gsub("vcf","tab",input) }
out<-sub(".gz","",out)
print(paste0("Output file name : ", out))

###################
#filter PASS
PASS<-unlist(strsplit(PASS,","))
print(paste("PASS=",PASS))
passvcf=gsub("vcf","pass.vcf",input)
passvcf=gsub("vcf$","vcf.gz",passvcf)
write.vcf( vcf[ vcf@fix[,'FILTER'] %in% PASS ], file=passvcf )

###################
#run annovar
annovarBin<-paste0( gsub("table_annovar.pl","",annovarBinPath), "/table_annovar.pl")
command=paste(annovarBin,"--buildver", reference , "--thread", threads, "--vcfinput --onetranscript --remove --otherinfo --protocol", protocols, "-operation", operations, passvcf, annovarDBpath , sep=" ")
print("Run table_annovar.pl")
print(command)
err<-system(command, ignore.stderr=T)
if(err>0){stop(err)}

################
#Functions

mergeAnnovarFiles<-function(dir="./"){
  
  variants<-list.files(path = dir)
  avinput<-list.files(path = dir, pattern="avinput")
  avinput<-fread(avinput[1])[,-c(6,7,8)]
  print(nrow(avinput))
  colnames(avinput)<-c(fixNames,colNames)
  nbfields<-ncol(avinput)
  
  for( variant in variants){
    dbname=sub(".*vcf\\.", "", variant)
    dbname=sub(".variant.*", "", dbname)
    dbname=sub(reference, "", dbname)
    dbname=sub("_","",dbname)
    dbname=sub(".txt", "", dbname)
    dbname=sub("gz.","",dbname)
    if( dbname %in% c("refGene.exonic","knownGene.exonic","ensGene.exonic") ){
      print(paste0("ANNOTATING : ",dbname))
      tmp<-fread(variant)[,1:8]
      colnames(tmp)<-c( paste0("GeneDetail.", dbname), paste0("ExonicFunc.", dbname), paste0("AAChange.", dbname), fixNames )
      avinput<-merge( tmp, avinput, by=fixNames, all.y=TRUE )
    }
    
    if( dbname %in% c("refGene","knownGene","ensGene") ){
      print(paste0("ANNOTATING : ",dbname))
      tmp<-fread(variant)[,1:7]
      colnames(tmp)<-c( paste0("Func.", dbname), paste0("Gene.",dbname), fixNames )
      avinput<-merge( tmp, avinput, by=fixNames, all.y=TRUE )
    }
    
    if( dbname %in% c("cytoBand","genomicSuperDups")){
      print(paste0("ANNOTATING : ",dbname))
      tmp<-fread(variant)[,2:7]
      colnames(tmp)<-c( dbname, fixNames)
      avinput<-merge( tmp, avinput, by=fixNames, all.y=TRUE )
    }
    
    if( dbname %in% c("multianno") ){
      print(paste0("ANNOTATING : ",dbname))
      suppressWarnings(tmp<-fread(variant))
      tmp<-tmp[,1:(ncol(tmp)-nbfields+2)]
      colnames(tmp)[1:5]<-fixNames
      avinput<-merge( tmp, avinput, by=fixNames, all.y=TRUE )
    }
    
  }
  return(avinput)
  
}

##########################################
# Retrieve Strand from refGene file
# It cannot work if refGene is not used !
getStrand<-function(avtmp){
  
  if ( is.null(avtmp$Gene.refGene)) { return(avtmp) }
  refGene<-fread( paste0(annovarDBpath, "/", reference, "_refGene.txt") )
  avtmp$symbol<-gsub("\\(.*","",avtmp$Gene.refGene)
  avtmp[ ,symbol2:=tstrsplit(symbol,";",keep=1) ]
  tmp<-unique(refGene[,c(13,4)])
  colnames(tmp)<-c("symbol","Strand")
  dup<-tmp$symbol[duplicated(tmp$symbol)]
  tmp[ symbol %in% dup,]$Strand<-"+/-"
  tmp<-tmp[ !duplicated(tmp$symbol), ]
  avtmp<-merge(avtmp,tmp, by="symbol", all.x=T)
  avtmp[,symbol:=NULL]
  avtmp[,symbol2:=NULL]
  return(avtmp)
}

##########################################
# Retrieve Context from fasta reference file
getContextAnnotation<-function(avtmp){
  
  reffile<-list.files( path = annovarDBpath, pattern=paste0(reference,".fa"), full.names =T )
  ref<-readDNAStringSet(reffile)
  avtmp$context=apply( avtmp, 1, function(x) getContext(ref,x["CHROM"],x["POS"],10) )
  avtmp$trinucleotide_context=apply( avtmp, 1, function(x) getContext(ref,x["CHROM"],x["POS"],1) )
  avtmp$trinucleotide_context<-sub("(.).(.)","\\1x\\2",avtmp$trinucleotide_context)  
  return(avtmp)
}

getContext<-function(ref,chr,pos,win){
  
  ctx<-subseq(ref[[chr]], max( 1, as.numeric(pos)-win ), min( length(ref[[chr]]), as.numeric(pos)+win ))
  return(as.character(ctx))
}


##################
#PostProcess
print("Merge annovar outputs")
avoutput<-mergeAnnovarFiles()
print(nrow(avoutput))
print("Strand annotation")
avoutput<-getStrand(avoutput)
print(nrow(avoutput))
print("Context annotation")
avoutput<-getContextAnnotation(avoutput)
print(nrow(avoutput))
print("Reorder columns")
dbnames<-colnames(avoutput)[ ! colnames(avoutput) %in% c(fixNames,"Strand","context","trinucleotide_context",colNames) ]
setcolorder(avoutput, c(fixNames,dbnames,"Strand","context","trinucleotide_context",colNames) )
print(nrow(avoutput))
print("Order by chromosomic location")
avoutput$num=avoutput$Chr
avoutput[ num=="chrX", num :="23"]
avoutput[ num=="chrY", num :="24"]
avoutput[, num:=as.integer(gsub("chr","",num))]
avoutput<-avoutput[ order(c(num,Start,End,Ref,Alt)) ]
avoutput[, num:=NULL]
print("Write output")
write.table(avoutput, file=out, row.names=F, sep="\t", quote=F)
print("The end !")
print(nrow(avoutput))
 
  

