#!/usr/bin/env Rscript

## Takes the final_tus files from any number of inputs and converts it 
## into consensus calls while reporting some statistics

## Requires a tuFilePath file which contains the paths to the final TU calls in column 1 and the experimental 
## design in columns 2..N. Label columns with whatever labels you choose.

## Read in arguments
suppressMessages(library(argparse))
parser <- ArgumentParser(description="Select transcriptional models based on data.")
parser$add_argument("pathList",help="File containing all paths to final_tus.txt.")
parser$add_argument("out",help="Directory to write consensus calls to") 

args <- parser$parse_args()

source("consensusLib.R")
dir.create(args$out)
dir.create(file.path(args$out,"qc_plots"))
pFile=fread(args$pathList)

## Save column names and rename
oldNames=colnames(pFile)
pFile = renameCols(pFile)

## Determine which rows are replicates
pFile=findRepl(pFile)

## Assign unique ids to each row just to protect against later rearrangement
pFile[,uid:=1:nrow(pFile)]

## Read in all files
tus=list()
for(i in sort(unique(pFile$uid))){
    tus[[i]]=fread(pFile[uid==i]$path,colClasses=c("character","numeric","numeric","character","character","factor","character","character","numeric","numeric","numeric","character","character"))
}

##
## TUs
##

## Look at concordance between replicates, including untranslated calls
## Note: Genes that are inactive in both replicates are not counted
concord.rpl=list()
for(i in sort(unique(pFile$rid))){
    tus.rpl=tus[pFile$uid[pFile$rid==i]]
    concord.rpl[[i]] <- calcTuConcordance(tus.rpl,include.untrans=TRUE)
}

## Roughly 84% of all calls are concordant between replicates including untranslated vs translated
concord.rpl.dist=table(rbindlist(concord.rpl)$per)
concord.rpl.dist/sum(concord.rpl.dist)
mean.tu.rpl.agree=rbindlist(concord.rpl)[,mean(per),by=GENEID]
setnames(mean.tu.rpl.agree,"V1","perTuRplAgree")
## concord.rpl.dist[1]/sum(concord.rpl.dist)

## Look at concordance between all runs, ignoring untranslated genes 
concord.all <- calcTuConcordance(tus)
## hist(concord.all$per)
setnames(concord.all,"per","perTuAgree")


##
## Runthrough
##

## Look at concordance between replicates for runthrough, 
## Note: Genes that are inactive in both replicates are not counted
concord.rpl.rt=list()
for(i in sort(unique(pFile$rid))){
    tus.rpl.rt=tus[pFile$uid[pFile$rid==i]]
    concord.rpl.rt[[i]] <- calcRtConcordance(tus.rpl.rt)
}
concord.rpl.dist.rt=table(rbindlist(concord.rpl.rt)$per)
concord.rpl.dist.rt/sum(concord.rpl.dist.rt)
ag.rt=rbindlist(concord.rpl.rt)[,agree:=(per==0 | per==1),]
mean.rt.rpl.agree=ag.rt[,mean(agree),by=GENEID]
setnames(mean.rt.rpl.agree,"V1","rtRplAgree")


## Look at concordance across all runs
concord.all.rt <- calcRtConcordance(tus)
## hist(concord.all.rt$per,breaks=0:length(tus)/length(tus))
setnames(concord.all.rt,"per","percent.rt")

##
## 5' and 3' ends
##

concord.rpl.co=list()
for(i in sort(unique(pFile$rid))){
    tus.rpl.co=tus[pFile$uid[pFile$rid==i]]
    concord.rpl.co[[i]] <- calcCoordConcordance(tus.rpl.rt)
}
concord.rpl.all=rbindlist(concord.rpl.co)
##hist(concord.rpl.all$fMAD,breaks=100)
## Fraction of Tus that have different starts between replicates
## sum(concord.rpl.all$fMAD==0,na.rm=TRUE)/sum(!is.na(concord.rpl.all$fMAD))
## sum(concord.rpl.all$tMAD==0,na.rm=TRUE)/sum(!is.na(concord.rpl.all$tMAD))
## Quantile distiribtion of MAD
## quantile(concord.rpl.all$fMAD,na.rm=TRUE,c(seq(0,0.9,0.05),seq(0.9,1,0.01)))
## quantile(concord.rpl.all$tMAD,na.rm=TRUE,c(seq(0,0.9,0.05),seq(0.9,1,0.01)))


# concord.rpl.dist.rt/sum(concord.rpl.dist.rt)

## Look at concordance across all runs
concord.all.co <- calcCoordConcordance(tus)

## Compute Mean GOF for each transcript
source("consensusLib.R")
gof.mean=calcPerTxMeanGof(tus)
setnames(gof.mean,"TXNAME","consensus")

##
## Merge all of these into one big table that can then be used (along with GOF) to pick consensus call
##
big=merge(concord.all,mean.tu.rpl.agree,by="GENEID",all=TRUE)
big=merge(big,concord.all.rt,by="GENEID",all=TRUE)
big=merge(big,mean.rt.rpl.agree,by="GENEID",all=TRUE) 
big=merge(big,concord.all.co,by="GENEID",all.x=TRUE)
big=merge(big,gof.mean,by=c("GENEID","consensus"),all.x=TRUE)
big=merge(big,rbindlist(tus)[,c("GENEID","TXTYPE"),with=FALSE],by="GENEID",all.x=TRUE)

## Filter by some criteria
## Either the calls are highly concordant or they tend to have very close beginnings and endings
big.filter=big[((perTuAgree>=0.8 & perTuRplAgree>=0.5) | (fMAD<=500 & tMAD <=500)) & percent.rt<=0.5,]

## Print out some QC plots
qcPlots(big,file.path(args$out,"qc_plots"),"unfiltered")
qcPlots(big.filter,file.path(args$out,"qc_plots"),"filtered")


### Produce consensus calls
tx.db=unique(rbindlist(tus)[,1:8,with=FALSE])
con.tx=tx.db[TXNAME %in% big.filter$consensus]
setnames(con.tx,"TXNAME","TXID")
con.tx[,GENEID:=NULL]
setorder(con.tx,"TXCHROM","TXSTART","TXEND","GENENAME","TXID","TXTYPE","TXSTRAND")

## Write output
write.table(big,file.path(args$out,"consensus_stats.txt"),row.names=FALSE,quote=FALSE,sep="\t")
write.table(con.tx,file.path(args$out,"consensus_tus.txt"),row.names=FALSE,quote=FALSE,sep="\t")
