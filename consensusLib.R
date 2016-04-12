#!/usr/bin/env Rscript

## Library of functions to do consensus calling from multiple TuSelecter calls

args=list()
args$pathList="/sonas-hs/siepel/hpc_norepl/home/ndukler/Projects/celastrol/results/tu_selecter_new/finalTuPath.txt"
args$parallel=2
args$out="/sonas-hs/siepel/hpc_norepl/home/ndukler/Projects/celastrol/results/tu_selecter_new/consensus_tus"

library(data.table)
library(doParallel)
library(ggplot2)

## Rename columns so later code can assume column names
renameCols <- function(pFile){
    ## Rename columns for ease of use
    if(ncol(pFile)>1){
        setnames(pFile,oldNames,c("path",paste0("V",1:(length(oldNames)-1))))
    } else {
        setnames(pFile,oldNames,c("path"))
    }
    return(pFile)
}

## Determine which rows are replicates
findRepl <- function(pFile){
    ind=unique(pFile[,-1,with=FALSE])
    ind[,rid:=1:nrow(ind)]
    return(merge(pFile,ind,by=colnames(pFile[,-1,with=FALSE])))
}

## Wrapper function for merge
mymerge <- function(x,y) {
    return(merge(x,y,by="GENEID",suffixes=c(".x",".y")))
}

## Given a list of calls look at concordance
calcTuConcordance <- function(sub.tus,include.untrans=FALSE,cores=args$parallel){
    ## Bind together names for each replicate
    temp=lapply(sub.tus,function(x) return (x[,c("GENEID","TXNAME"),with=FALSE]))
    for(i in 1:length(temp)){
        setnames(temp[[i]],colnames(temp[[i]])[2],paste0("V",i))    
    }
    txnames=Reduce(mymerge,temp)
    ## Remove genes that are always untranscribed
    include=which(apply(txnames[,-1,with=FALSE]=="Untranscribed",1,sum)<(ncol(txnames)-1))
    ## Set up data.table to not count Untranscribed genes if that option is enabled
    txnames=txnames[include]
    if(include.untrans==FALSE){
        txnames[txnames=="Untranscribed"]=NA
    }

    ## Split data by core
    alloc=rep_len(x=1:cores, length.out=nrow(txnames))
    split.tx=list()
    for(i in 1:cores){
        split.tx[[i]]=txnames[alloc==i]
    }

    con=rbindlist(mclapply(split.tx,function(x){
        out=(rbindlist(apply(x,1,function(k) {
            z=table(k[-1],useNA="no")
            return(data.table(GENEID=k[1],consensus=names(z)[which.max(z)],per=max(z)/sum(z)))
        })))
        return(out)
      },mc.cores=cores))
    return(con)
}


## Computes how often RT calls agree based on a given alpha level
calcRtConcordance <- function(sub.tus,alpha=0.05){
    ## Bind together names for each replicate
    temp=lapply(sub.tus,function(x) return (x[,c("GENEID","q"),with=FALSE]))
    for(i in 1:length(temp)){
        ## Convert all catagorical variables to 1="RT" or 0="non-RT"
        temp[[i]]$q[temp[[i]]$q=="non-RT"]=0
        temp[[i]]$q[temp[[i]]$q=="RT"]=1
        temp[[i]]$q=as.numeric(temp[[i]]$q)
        setnames(temp[[i]],colnames(temp[[i]])[2],paste0("V",i))
    }
    rt=Reduce(mymerge,temp)
    ## Remove genes that are always untranscribed
    include=which(apply(is.na(rt[,-1,with=FALSE]),1,sum)<(ncol(rt)-1))
    ## Remove genes that are always untranscribed
    rt=rt[include]

    ## How many RT calls are there per row
    hold=rt$GENEID
    rt[rt<=alpha]=FALSE
    rt[rt>alpha]=TRUE
    rt$GENEID=hold
    out=data.table(GENEID=rt$GENEID,per=apply(rt[,-1,with=FALSE],1,sum,na.rm=TRUE))
    out[,per:=per/(ncol(rt)-1)]
    return(out)
}

## Computes mean absolute deviation for the 5' and 3' ends of genes
## Ignores untranscribed calls 
calcCoordConcordance <- function(sub.tus){
    ## Get gene ids for each row
    info=data.table(GENEID=sub.tus[[1]]$GENEID,TXSTRAND=sub.tus[[1]]$TXSTRAND)
    ## Bind together names for each replicate
    txname=as.data.table(do.call("cbind", lapply(sub.tus,function(x) return (x$TXNAME))))
    ## Bind together txstart for each replicate
    start=as.data.table(do.call("cbind", lapply(sub.tus,function(x) return (x$TXSTART))))
    ## Bind together txend ends for each replicate
    end=as.data.table(do.call("cbind", lapply(sub.tus,function(x) return (x$TXEND))))

    ## Set starts and ends of untranscribed transcripts to NA
    start[txname=="Untranscribed"]=NA
    end[txname=="Untranscribed"]=NA

    ## Compute variance in starts and ends for all transcripts
    sMAD=data.table(row=1:nrow(info),info,mad=apply(start,1,function(x) mean(abs(x-mean(x,na.rm=TRUE)),na.rm=TRUE)))
    eMAD=data.table(row=1:nrow(info),info,mad=apply(end,1,function(x) mean(abs(x-mean(x,na.rm=TRUE)),na.rm=TRUE)))
   
    ## Transform start & end to 3' and 5' with respect to strand
    fMAD=rbind(sMAD[TXSTRAND=="+"],eMAD[TXSTRAND=="-"])
    fMAD=fMAD[order(row)]
    fMAD[,row:=NULL]
    tMAD=rbind(sMAD[TXSTRAND=="-"],eMAD[TXSTRAND=="+"])
    tMAD=tMAD[order(row)]
    tMAD[,row:=NULL]

    return(data.table(GENEID=fMAD$GENEID,fMAD=fMAD$mad,tMAD=tMAD$mad))   
}

## Compute the mean reported GOF for each transcript id 
calcPerTxMeanGof <- function(tus){
    tx.all=rbindlist(tus)
    tx.gof=tx.all[,mean(gof),by=list(TXNAME,GENEID)]
    setnames(tx.gof,"V1","mean.gof")
    return(tx.gof)
}

## Produce some graphs and logs to give people an idea of the quality of their output
## in the unfiltered set.
qcPlots <- function(qc.table,out,label){
## First a graph of the average percenteage of tus that agree
    pdf(file.path(out,paste0("percentTuAgree_",label,".pdf")))
    ggplot(qc.table,aes(x=perTuAgree))+
        geom_histogram()
    dev.off

    pdf(file.path(out,paste0("percentTuAgreeByTxtype_",label,".pdf")))
    ggplot(qc.table,aes(x=perTuAgree))+
        facet_wrap(~TXTYPE)+
            geom_histogram()
    dev.off()
    
    ## Then show mean GOF
    pdf(file.path(out,paste0("gof.pdf_",label,".pdf")))
    ggplot(qc.table,aes(x=mean.gof))+
    geom_histogram()
    dev.off()
    
    pdf(file.path(out,paste0("gofByTxtype_",label,".pdf")))
    ggplot(qc.table,aes(x=mean.gof))+
        facet_wrap(~TXTYPE)+
            geom_histogram()
    dev.off()
    
    ## Show the percent with RT
    pdf(file.path(out,paste0("runthroughtByTxtype_",label,".pdf")))
    ggplot(qc.table,aes(x=percent.rt))+
    geom_histogram()
    dev.off()
    
    pdf(file.path(out,paste0("runthroughByTxtype_",label,".pdf")))
    ggplot(qc.table,aes(x=percent.rt))+
        facet_wrap(~TXTYPE)+
            geom_histogram()
    dev.off()    
}
