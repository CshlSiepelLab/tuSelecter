#!/usr/bin/env Rscript

## Reimplementation of TU selecter, hopefully much more efficient and gets rid of bugs
## Read in arguments
suppressMessages(library(argparse))
parser <- ArgumentParser(description="Select transcriptional models based on data.")
parser$add_argument("txtable",help="Path to the table of transcriptional units.") 
parser$add_argument("bwp", help="Bigwig file of proseq data for the plus strand.");
parser$add_argument("bwm", help="Bigwig file of proseq data for the minus strand.");
parser$add_argument("-w","--up-window",type='integer',help="How close an upstream active gene must be to check for runthrough transcription (default 10^5 bp).",default=10^5)
parser$add_argument("-t","--tile",type='integer',help="Length of tile used for tiling genome.(Default 100bp)",default=100)
parser$add_argument('-p',"--parallel",type='integer',help="Number of cores to use",default=1)
parser$add_argument("-o","--out",type='character',help="Where output from the analysis goes",default=".")

args <- parser$parse_args()

log=file.path(args$out,"log.txt")
args$min_dist=3*args$tile

source("rtHMMLib.R")

## Set-up directory structure
dir.create(args$out)

## Read in transcript models
write("Reading in transcript models...",file=log)
txTable=readTxTable(args$txtable)
## Filter out transcripts that are too small to be analyzed
write("Filtering transcript models...",file=log,append=TRUE)
txTable.filtered=filterTxTable(txTable,args$tile)
txTable.filtered.gr=txTableToGR(txTable.filtered)

## Read in bigWigs
write("Reading in bigwigs...",file=log,append=TRUE)
bw=readBw(args$bwp,args$bwm)

###############
## TU SELECTION
###############

## Construct bounded region covered by all TxModels per gene
geneTab=boundTxModels(txTable.filtered)

## Convert geneTab models into bins
write("Tiling transcript models...",file=log,append=TRUE)
geneTiles=tileTus(geneTab,tile=args$tile)
geneTiles.gr=tileToGR(geneTiles)
rm(list=c("geneTiles"))

## Get sums of reads in bins
write("Summing reads over tiled transcript models...",file=log,append=TRUE)
geneTileSums=sumBwOverGR(geneTiles.gr,bw)

## Cleanup gene tile
rm(list=c("geneTiles.gr"))

## Initialize parallel environment
## initParallel(args$parallel)

## Set up parameters for TU selection
p.grid=seq(0.3,1,.05)
p.grid.back=seq(0.01,.03,.01)
trans.params=expand.grid(p.grid)
back.params=expand.grid(p.grid.back)

uID=unique(txTable.filtered$ensg.id)
## Split uID into chunks
split.uID=split(uID, ceiling(seq_along(uID)/(length(uID)/args$parallel)))

## Run TU selection in parallel
write("Selecting TUs (this may take a while)...",file=log,append=TRUE)
system.time(tuList <- rbindlist(mclapply(split.uID, function(x){
    out=rbindlist(apply(as.matrix(x),1,function(x) geneTuSelect(ensg.id=x,txTable.gr=txTable.filtered.gr,
        geneTileSums,trans.params,back.params)))
},mc.preschedule=FALSE, mc.cores = args$parallel)))

## Write out intermediate file
write.table(tuList,file.path(args$out,"tuList.txt"),row.names=FALSE,quote=FALSE)

pdf(file.path(args$out,"gof_stats.pdf"))
    plot(tuList$gof[tuList$gof!=-1],tuList$lnmr[tuList$gof!=-1],main="Goodness of fit statistics",xlab="GOF",ylab="Longest Non-matching Run")
dev.off()

## Plot some GOF statistics
pdf(file.path(args$out,"gof_stats_cdf_transcribed.pdf"))
    plot(ecdf(tuList$gof[tuList$gof>0 & !is.na(tuList$TXCHROM)]),main="CDF GOF (transcribed genes only)",xlab="GOF")
    plot(ecdf(tuList$lnmr[tuList$gof>0 & !is.na(tuList$TXCHROM)]),main="CDF LNMR (transcribed genes only)",xlab="LNMR")
dev.off()

pdf(file.path(args$out,"gof_cdf_all.pdf"))
    plot(ecdf(tuList$gof[tuList$gof>0]),main="CDF GOF (all genes)",xlab="GOF")
    plot(ecdf(tuList$lnmr[tuList$gof>0]),main="CDF LNMR (all genes)",xlab="LNMR")
dev.off()

## Release memory for variables no longer in use
rm(list=c("split.uID","uID","geneTileSums"))

#######################
## RUNTHROUGH DETECTION
#######################
write("Beginning runthrough detection algorithm...",file=log,append=TRUE)
nrow(tuList)

source("rtHMMLib.R")
## Split transcripts by strand and chromosome
ssTu=splitSortTuList(tuList)

sum(unlist(lapply(unlist(ssTu,FALSE),nrow)))==nrow(tuList)

## Sanity checks
## sum(duplicated(unlist(lapply(unlist(ssTu,FALSE),function(x) return(x$GENEID)))))
## sum(duplicated(unlist(lapply(unlist(ssTu,FALSE),function(x) return(x$LINK)))))
## sum(unlist(lapply(unlist(ssTu,FALSE),function(x) return(x$GENEID!=x$LINK))))

## Find NAT for all transcripts
write("Finding nearest active transcripts for each transcript...",file=log,append=TRUE)
nat.p <- mclapply(ssTu$p, function(x) findNAT(x),mc.cores=args$parallel)
nat.m <- mclapply(ssTu$m, function(x) findNAT(x),mc.cores=args$parallel)

## Sanity checks
## for(i in 1:length(ssTu$p)){print(sum(ssTu$p[[i]]$TXNAME=="Untranscribed")+nrow(nat.p[[i]])==nrow(ssTu$p[[i]]))}
## for(i in 1:length(ssTu$m)){print(sum(ssTu$m[[i]]$TXNAME=="Untranscribed")+nrow(nat.m[[i]])==nrow(ssTu$m[[i]]))}
## sum(duplicated(unlist(nat.p,FALSE)$LINK))
## sum(duplicated(unlist(nat.m,FALSE)$LINK))
## sum(unlist(lapply(nat.p,nrow)))+sum(unlist(lapply(nat.m,nrow)))+sum(tuList$TXNAME=="Untranscribed")
## nrow(tuList)

## Cleanup ssTu
rm(list=c("ssTu"))

## Convert NATs and TUs to the respective RT regions
## Where dist(NAT,TUs)<=minDist, TU is added to auto-RT list
rtTUs=natToTu(nat.p,nat.m,tuList)

## Cleanup NATs
rm(list=c("nat.p","nat.m"))

## Tile RT transcripts
write("Tiling RT transcript models...",file=log,append=TRUE)
rt.tile=tileTus(rbind(rtTUs$p,rtTUs$m),tile=50)
rt.tile.gr=tileToGR(rt.tile)
rm(list=c("rt.tile"))

## Get reads in tiles
write("Summing reads over tiles...",file=log,append=TRUE)
rtTileSums=sumBwOverGR(rt.tile.gr,bw)

## Cleanup RT.tile.gr
rm(list=c("rt.tile.gr"))

rt.split=splitRtTask(list(rtTileSums$p,rtTileSums$m))

## Do garbage collection to reduce memory use when forking
gc()

## Compute runthrough
print("Computing runthough scores...",file=log,append=TRUE)
rt.scores <- parallel.rt(rt.split)

## Correct for multiple testing
rt.scores[,q:=p.adjust(rt.p)]
## Reincorporate the tus that were automatically labeled as runthrough or not runthrough
rt.all=rbind(rt.scores,data.table(GENEID=rtTUs$auto$GENEID,rt.p="RT",q="RT"),data.table(GENEID=rtTUs$auto.non$GENEID,rt.p="non-RT",q="non-RT"))

## Merge RT calls with TU_selections
tuList.all=merge(tuList,rt.all,by="GENEID",all.x=TRUE)
## Set column order, do some cleanup and sorting
setcolorder(tuList.all,c(colnames(tuList),"rt.p","q"))
tuList.all[,LINK:=NULL]
tuList.all=tuList.all[order(TXCHROM,TXSTART,TXEND)]

## Write out all TUs and their RT scores
write.table(tuList.all,file=file.path(args$out,"final_tus.txt"),row.names=FALSE,quote=FALSE)

write("Done!",file=log,append=TRUE)
