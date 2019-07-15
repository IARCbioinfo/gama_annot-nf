#!/usr/bin/env Rscript
library(data.table)

somatic<-list.files(pattern="*somatic.snvs.fix.tab$")
germline<-list.files(pattern="*germline.fix.tab$")
mutect2<-list.files(pattern="*calls.tab$")

#strelka germline
print("strelka germline")
for (file in germline){

#TODO : rewrite this part
	print(file)

	vcf<-fread(file)
        vcf<-vcf[vcf$FILTER=="PASS" & ALT!="."]

	T<-colnames(vcf)[10] # name of sample

	if (nrow(vcf)>0){
		vcf[ ALT!="."  ,  COV :=  tstrsplit( tstrsplit(get(T),":")[4][[1]], "," )  ]

		#snvs
		vcf[ grep("SB", vcf$FORMAT)  ,  VAF_REF:=  tstrsplit( tstrsplit(get(T),":")[6][[1]], "," )[1]  ]
		vcf[ grep("SB", vcf$FORMAT)  ,  VAF_ALT:=  tstrsplit( tstrsplit(get(T),":")[6][[1]], "," )[2]  ]

		#indel
		vcf[ grep("DPI", vcf$FORMAT), VAF_REF:=  tstrsplit( tstrsplit(get(T),":")[5][[1]], "," )[1] ]
		vcf[ grep("DPI", vcf$FORMAT), VAF_ALT:=  tstrsplit( tstrsplit(get(T),":")[5][[1]], "," )[2] ]

		vcf$VAF_REF<-as.numeric(vcf$VAF_REF)/as.numeric(vcf$COV)
		vcf$VAF_ALT<-as.numeric(vcf$VAF_ALT)/as.numeric(vcf$COV)
	}

        filename<-gsub("_L00.","",file)
        filename<-gsub("BQSrecalibrated.bam","",filename)
	filename<-gsub(".fix",".pass",filename)
        print(filename)
	write.table(vcf, file=filename , sep="\t", row.names=F, quote=F)
}


#strelka somatic
print("strelka somatic")
for (file in somatic){
	print(file)
	
	snvs<-fread(file)
	snvs<-snvs[snvs$FILTER=="PASS" & ALT!="."]

	dict<-list( "A"=6, "C"=7, "G"=8, "T"=9 )

	#snvs
	if(nrow(snvs)>0){

		snvs[  ,  index:=dict[ ALT ] ]
		snvs$index<-as.integer(snvs$index)
		
		snvs[  ,  VAF_N  :=sapply( 1:nrow(snvs), function(x) { sub(",.*","",strsplit(snvs$NORMAL[x],":")[[1]][snvs$index[x]]) } ) ]
		snvs[  ,  VAF_T  :=sapply( 1:nrow(snvs), function(x) { sub(",.*","",strsplit(snvs$TUMOR[x],":")[[1]][snvs$index[x]]) } ) ]
		
		snvs[ , Cov_N :=  tstrsplit(NORMAL,":",keep=2)  ]
		snvs[ , Cov_T :=  tstrsplit(TUMOR,":",keep=2)  ]

		snvs[  ,  Cov_alt_N  :=sapply( 1:nrow(snvs), function(x) { sub(",.*","",strsplit(snvs$NORMAL[x],":")[[1]][snvs$index[x]]) } ) ]
		snvs[  ,  Cov_alt_T  :=sapply( 1:nrow(snvs), function(x) { sub(",.*","",strsplit(snvs$TUMOR[x],":")[[1]][snvs$index[x]]) } ) ]

		snvs$VAF_N<-as.numeric(snvs$VAF_N)/as.numeric(snvs$Cov_N)
                snvs$VAF_T<-as.numeric(snvs$VAF_T)/as.numeric(snvs$Cov_T)

                snvs<-snvs[ , index:=NULL ]
	}

	#indels
	print(sub("snvs","indels",file))
	indels<-fread(sub("snvs","indels",file))
	indels<-indels[indels$FILTER=="PASS" & ALT!="."]

	if(nrow(indels)>0){
		

		indels[ , VAF_N := tstrsplit(NORMAL,":",keep=5) ]
		indels[ , VAF_T := tstrsplit(TUMOR,":",keep=5) ]
		indels[ , VAF_N := sub(",.*","",VAF_N) ]
		indels[ , VAF_T := sub(",.*","",VAF_T) ]
		
		indels[ , Cov_N :=  tstrsplit(NORMAL,":",keep=2)  ]
		indels[ , Cov_T :=  tstrsplit(TUMOR,":",keep=2)  ]

		indels[ , Cov_alt_N := tstrsplit(NORMAL,":",keep=5) ]
		indels[ , Cov_alt_T := tstrsplit(TUMOR,":",keep=5) ]
		indels[ , Cov_alt_N := sub(",.*","",Cov_alt_N) ]
		indels[ , Cov_alt_T := sub(",.*","",Cov_alt_T) ]

		indels$VAF_N<-as.numeric(indels$VAF_N)/as.numeric(indels$Cov_N)
                indels$VAF_T<-as.numeric(indels$VAF_T)/as.numeric(indels$Cov_T)

	}

        tab <-snvs
	if ( ncol(snvs) == ncol(indels) ){ tab<-rbind(snvs,indels) }

	filename<-gsub("_L00._","",file)
        filename<-gsub("BQSrecalibrated.bam","",filename)
	filename<-gsub("snvs.","",filename)
	write.table(tab, file=filename , sep="\t", row.names=F, quote=F)
}


#Mutect2
print("mutect2 somatic")
for ( file in mutect2 ){

        print(file)

        snvs<-fread(file)

        Cov_N<-vcf[ , tstrsplit( NORMAL, ":", keep = 2 ) ][ , tstrsplit( V1, "," ) ]
        Cov_alt_N<-Cov_N[ , as.numeric(V2) ]
        Cov_N<-Cov_N[ , as.numeric(V1)+as.numeric(V2) ]

        VAF_N<- vcf[ , c("VAF_N") := tstrsplit( NORMAL, ":", keep = 3 ) ][,VAF_N]

        Cov_T<-vcf[ , tstrsplit( TUMOR, ":", keep = 2 ) ][ ,  tstrsplit( V1, "," ) ]
        Cov_alt_T<-Cov_T[ , as.numeric(V2) ]
        Cov_T<-Cov_T[ , as.numeric(V1)+as.numeric(V2) ]

        VAF_T<- vcf[ , c("VAF_T") := tstrsplit( TUMOR, ":", keep = 3 ) ][,VAF_T]

        foo<-cbind(vcf, Cov_N, Cov_T, VAF_N, VAF_T, Cov_alt_N, Cov_alt_T)
        
        filename<-gsub("_L00._","",file)
        filename<-gsub("_BQSrecalibrated","",filename)
        filename<-gsub("_calls", "_pass",filename)
        write.table(foo[FILTER=="PASS",], file=filename, sep="\t", row.names=F, quote=F)
}
