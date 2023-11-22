#! /usr/bin/env Rscript
#####################################################################################
#
# Title		: gatContext.r
# Author	: CahaisV@iarc.fr
# Date		: 31/07/2020
# Last Update	: 31/07/2020
#
#####################################################################################

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GetoptLong))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(parallel))

userAnnot=""
annovarDBpath="/data/databases/annovar/mm10db"

GetoptLong(matrix(c("annovarDBpath|a=s",  "path to annovarDB",
                    "userAnnot|u=s",      "additional user annotations"
), ncol=2, byrow=TRUE))

avinput<-list.files(pattern="tsv$")
if(!file.exists(avinput)){ stop("ERROR : Annovar output not found") }
if(file.info(avinput)$size<100){ stop("ERROR : Annovar output empty") }
out=gsub(".tsv",".1.tsv",avinput[1])
print(paste0("Output file name : ", out))
avinput<-fread(avinput)
startcols<-colnames(avinput)[1:(which(colnames(avinput)=="CHROM")-1)]
endcols<-colnames(avinput)[which(colnames(avinput)=="CHROM"):length(avinput)]

annovarDB=basename(annovarDBpath)
reference=gsub("db","",annovarDB)

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
  #avtmp$context=apply( avtmp, 1, function(x) getContext(ref,x["CHROM"],x["POS"],10) )
  #avtmp$trinucleotide_context=apply( avtmp, 1, function(x) getContext(ref,x["CHROM"],x["POS"],1) )
  avtmp[, context:=getContext(ref,CHROM,POS,10)]
  avtmp$context<-as.character(avtmp$context)
  avtmp[, trinucleotide_context:=substr(context, 10, 12)]
  avtmp$trinucleotide_context<-sub("(.).(.)","\\1x\\2",avtmp$trinucleotide_context)  
  return(avtmp)
}

getContext<-function(ref,chr,pos,win){
  return ( mclapply(1:length(chr), function(x) toString( subseq(ref[[ chr[x] ]], max( 1, as.numeric(pos[x])-win ), min( length(ref[[ chr[x] ]]), as.numeric(pos[x])+win )) )) )
}

addUserAnnot<-function(avtmp,dt_reg){
  
  dt_pos <- avtmp[ , .( Chr, Start, End ) ]
  dt_pos <- unique( dt_pos )
  dt_pos[ , Start := Start - 1 ]
  setkey( dt_pos, Chr, Start, End )
  ## overlap
  dt_ovlap <- foverlaps( dt_reg, dt_pos )[ ! is.na( Start ) ]
  dt_ovlap <- unique( dt_ovlap[ , .( Chr, Start, End, type ) ] )
  dt_ovlap[ , Start := Start + 1 ]
  dt_ovlap <- dcast( dt_ovlap, Chr + Start + End ~ type, fill = 0, fun.aggregate = length, value.var = "type" )
  
  avtmp<-merge(avtmp,dt_ovlap,by=c("Chr","Start","End"), all.x=T)
  avtmp[is.na(avtmp)]<-0
  return(avtmp)
}

##################
#main

print("Strand annotation")
avoutput<-getStrand(avinput)
print(nrow(avoutput))

print("Context annotation")
avoutput<-getContextAnnotation(avoutput)
print(nrow(avoutput))

print("User annotation")
if (grepl("RData",userAnnot)){
    load(userAnnot)
    if (exists("dt_reg")){
      if (sum( colnames(dt_reg)==c("Chr","Start","End","type") )== 4){
        avoutput<-addUserAnnot(avoutput,dt_reg)
      } else {print(paste0("WARNING : column names of ", userAnnot, " must be Chr, Start, End and type. Cancel user annot"))}
    } else { print(paste0("WARNING : ", userAnnot, " must contain a dt_reg object. Cancel user annot"))}
} else if (userAnnot!=""){
    dt_reg<-fread(userAnnot)
    if (sum( colnames(dt_reg)==c("Chr","Start","End","type") )== 4){
        avoutput<-addUserAnnot(avoutput,dt_reg)
    }else {print(paste0("WARNING : column names of ", userAnnot, " must be Chr, Start, End and type. Cancel user annot"))}
}

#reorder lines
print("Order by chromosomic location")
avoutput$num=avoutput$Chr
avoutput[ num=="chrX", num :="23"]
avoutput[ num=="chrY", num :="24"]
avoutput[, num:=as.integer(gsub("chr","",num))]
avoutput<-avoutput[ order(cnum,Start,End,Ref,Alt), ]
avoutput[, num:=NULL]

#reorder columns
setcolorder(avoutput, c(startcols,"Strand","context","trinucleotide_context",endcols) )

print("Write output")
write.table(avoutput, file=out, row.names=F, sep="\t", quote=F)
print("The end !")
print(nrow(avoutput))
