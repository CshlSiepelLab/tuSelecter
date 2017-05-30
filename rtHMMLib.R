#!/usr/bin/env Rscript

## Import libraries
library(data.table)
library(depmixS4)
library(rtracklayer)
## library(Gviz)
library(doParallel)

## Initialize parallel cluster
initParallel <-  function(cores=NA){
    if(!is.na(cores)){
        cl = makeCluster(cores)
        registerDoParallel(cl)
    } else {
        cl = makeCluster(1)
        registerDoParallel(cl)
    }
    return(cl)
}

## Read in txtable and return it as a GRanges object
readTxTable <- function(filename){
    tuTable=fread(filename)
    coln=c("chr","start","end","ensg.id","enst.id","strand","hngc.id","txtype")
    if(ncol(tuTable)!=8){
        stop("TxTable incorrect number of columns. Columns must be ", paste(coln,collapse=", "),".")
    }
    setnames(tuTable,colnames(tuTable),coln)
    return(tuTable)
}

## Removes transcripts that do not meet requirement to be analyzed
filterTxTable <- function(txTable,tile,chrom){
    if(nrow(txTable[!(end-start>tile & chr %in% chrom)])>0){
        write("Filtering out un-modelable transcripts...",file=log,append=TRUE)
        write(paste(txTable[!(end-start>tile & chr %in% chrom)]$ensg.id,collapse="\n"),file=log,append=TRUE)
    }
    txTable=txTable[end-start>tile & chr %in% chrom]
    return(txTable)
}

txTableToGR <- function(txTable){
    return(with(txTable,GRanges(chr, IRanges(start=start,end=end,name=enst.id), strand=strand,ensg.id,txtype,hngc.id)))
}

## Takes in TxTable and gets lowest and highest coordinate for each gene
boundTxModels <- function(txTable){
    geneTab=txTable[,list(chr[1],min(start),max(end),strand[1]),by=ensg.id]
    setnames(geneTab,colnames(geneTab)[2:5],c("chr","start","end","strand"))
    return(geneTab)
}

## Read in bigWigs
readBw <- function(plus,minus){
    out=list()
    out[["+"]]=import(plus,as='RleList')
    out[["-"]]=import(minus,as='RleList')
    return(out)
}

## Gets ranges of bigwig files
getChromInfo <- function(bwPaths,which.chrom=NULL){
    rl <- BigWigSelection(IRangesList(IRanges(-1,-1)))
    chromInfo = lapply(import(con=bwPaths[1],selection=rl,as='RleList'),function(x){return(x@lengths)})
    ## If a specific set of chromosomes is specified only get info for that list of chromosomes
    if(!is.null(which.chrom)){
        chromInfo=chromInfo[names(chromInfo) %in% which.chrom]
    }
    if(length(bwPaths)>1){
        for(i in 2:length(bwPaths)){
            cur=lapply(import(con=as.character(bwPaths[i]),selection=rl,as='RleList'),function(x){return(x@lengths)})
            if(!is.null(which.chrom)){
                cur=cur[names(cur) %in% which.chrom]
            }
            if(sum(!unlist(Map("==",chromInfo,cur)))>0){
                stop(paste("Bigwig", bwPaths[i],"was built with a different chromInfo file than",  bwPaths[1]))
            }
        }
    }
    return(chromInfo)
}

## Takes in 5 column table, table must contain chr,start,end,strand,symbol. The symbol column can be in
## any position and have any name. It will be renamed as symbol however in the output. 
tileTus <- function(geneTab,tile){
    ## write("Tiling genes...",stdout())
    ## Figure out which column contains the symbol
    symCol=which(!colnames(geneTab) %in% c("chr","start","end","strand"))
    
    ## Split transcripts by strand
    geneTabP=geneTab[geneTab$strand=="+",]
    geneTabM=geneTab[geneTab$strand=="-",]

    ## Generate tiles on plus strand
    tileP=geneTabP[,list(chr,seq(start,end-tile,tile)),by=eval(colnames(geneTabP)[symCol])]
    setnames(tileP,colnames(tileP)[3],"start")
    tileP[,end:=start+eval(tile-1),]
    setnames(tileP,colnames(tileP)[1],"symbol")

    ## Generate tiles on minus strand
    tileM=geneTabM[,list(chr,rev(seq(end-tile,start,-tile))),by=eval(colnames(geneTabM)[symCol])]
    setnames(tileM,colnames(tileM)[3],"start")
    tileM[,end:=start+eval(tile-1),]
    setnames(tileM,colnames(tileM)[1],"symbol")
    
    return(list(plus=tileP, minus=tileM))
}

## Converts tiles to a GRanges object
tileToGR <- function(tus){
    ## Convert to granges objects
    tus.plus.gr=with(tus$plus, GRanges(chr, IRanges(start=start,end=end,name=symbol), strand="+"))
    tus.minus.gr=with(tus$minus, GRanges(chr, IRanges(start=start,end=end,name=symbol), strand="-"))       
    return(list(plus=tus.plus.gr,minus=tus.minus.gr))
}

## Sums reads in GRranges object
sumBwOverGR <- function(tu.list,bw.list){
    ## write("Summing reads...",stdout())
    ## Loop over chromosomes on both strands
    tu.list$plus$sum=-1
    tu.list$minus$sum=-1
    for(chr in unique(seqnames(tu.list$plus))) {
        tu.list$plus[seqnames(tu.list$plus) == chr]$sum=as.numeric(sum(Views(bw.list[["+"]][[as.character(chr)]], ranges(tu.list$plus[seqnames(tu.list$plus) == chr]))))      
    }
    for(chr in unique(seqnames(tu.list$minus))) {
        tu.list$minus[seqnames(tu.list$minus) == chr]$sum=abs(as.numeric(sum(Views(bw.list[["-"]][[as.character(chr)]], ranges(tu.list$minus[seqnames(tu.list$minus) == chr])))))                
    }
    return(tu.list)
}

## Computes an indicator function for each potential transcript relative to the gene model
## for each bin as transcribed or untranscribed.
calc.indicator <- function(models.gr,bins,bin.size){
    ## Find bins that overlap each transcript model by <= 50%
    overlaps = findOverlaps(models.gr,bins,minoverlap = 0.5*bin.size)
    ## Compute indicator function for bins
    indicator=list()
    for(i in unique(overlaps@from)){
        indicator[[i]]= list(trans=overlaps@to[overlaps@from==i],
                     un.trans=setdiff(unique(overlaps@to),overlaps@to[overlaps@from==i]))
    }
    return(indicator)
}

## Produces a likihood for all models
liklihood <- function(trans.params,back.params,mod,binned.data){
    ll.mat=liklihood.grid(trans.params,back.params,mod,binned.data)
    for(j in 1:length(ll.mat)){
        if(!is.matrix(ll.mat[[j]])){
            ll.mat[[j]]=matrix(ll.mat[[j]],ncol=length(ll.mat[[j]]),byrow=TRUE)
        }
    }
    ll=list(background=sum(apply(ll.mat$background,1,function(x) max(x)+log(sum(exp(x-max(x)))))),
        transcribed=sum(apply(ll.mat$transcribed,1,function(x) max(x)+log(sum(exp(x-max(x)))))))
    return(ll) 
}

## Returns the likilhood off all states at all positions
liklihood.grid <- function(trans.params,back.params,mod,bins){
    return(binom.ll(trans.params,back.params,mod,bins))
}

## Uses a binomial model to find the log liklihood
binom.ll <- function(trans.params,back.params,indic,bins){
if(length(bins[indic$trans])>0){
        transcribed=apply(as.matrix(trans.params),1,function(x){
            val=bins[indic$trans]
            a=dbinom(val>0,size=1,prob=x[1],log=TRUE)
            return(a-log(nrow(trans.params)))
        })
      } else{
          transcribed=matrix(-log(nrow(trans.params)),ncol=nrow(trans.params),nrow=1)
      }
      if(length(bins[indic$un.trans])>0){
          background=apply(as.matrix(back.params),1,function(x){
              val=bins[indic$un.trans]
              a=dbinom(val>0,size=1,prob=x[1],log=TRUE)
              return(a-log(nrow(back.params)))
          })
      } else{
          background=matrix(-log(nrow(back.params)),nrow=1,ncol=nrow(back.params))
      }
      return(list(background=background,transcribed=transcribed))
}

## Takes a set of indicator function and merges non-identifiable models and selects one
## enst.id (will preferentially select the protein coding one.
reduceModels <- function(txMods,txInd){
    uTx=unique(txInd)
    if(length(uTx)<length(txInd)){
        keep=numeric(length(uTx))
        for(i in 1:length(uTx)){
            same=which(unlist(lapply(txInd,function(x) identical(uTx[[i]],x))))
            if(sum(txMods$txtype[same]=="protein_coding")>=1){
                same=same[which(txMods$txtype[same]=="protein_coding")]
            }
            keep[i]=same[1]
        }
        return(list(txMods=txMods[keep,],txInd=txInd[keep]))
    } else {
        return(list(txMods=txMods,txInd=txInd))
    }
}

## Takes in a set of models and computes a posterior
compute.model.posterior <-  function(trans.params,back.params,indicator,binned.data){
    ## integrate out parameters
    posterior=list()
    for(z in 1:length(indicator)){
        if(names(indicator)[z]=="" || is.null(names(indicator)[z])){
            posterior[[z]]=liklihood(trans.params,back.params,indicator[[z]],binned.data)
        }
        else {
            posterior[[names(indicator)[z]]]=liklihood(trans.params,back.params,indicator[[z]],binned.data)
        }
    }
     
    ## Compute normalization factor
    norm.factor.trans=lapply(posterior,function(x) max(x$transcribed)+log(sum(exp(x$transcribed-max(x$transcribed)))))
    norm.factor.back=lapply(posterior,function(x) max(x$background)+log(sum(exp(x$background-max(x$background)))))
    
    ## Compute model evidence under uniform prior
    model.evid.trans = unlist(norm.factor.trans)-log(length(trans.params))
    model.evid.back = unlist(norm.factor.back)-log(length(back.params))
    model.evid=model.evid.trans+model.evid.back
    model.norm.factor= max(model.evid)+log(sum(exp(model.evid-max(model.evid))))
    model.posterior=exp(model.evid.trans+model.evid.back-model.norm.factor)
    return(model.posterior)
}

## Checks goodness of fit of the selected TU model
gof.map <- function(mod.map,val,indicator){
    indic=list(transcribed=c(1:length(val)),un.trans=c(1:length(val)))
    l.gof=liklihood.grid(trans.params,back.params,indic,val)        

    if(!is.matrix(l.gof$transcribed)){
        l.gof$transcribed=matrix(l.gof$transcribed,ncol=length(l.gof$transcribed),byrow=TRUE)
    }
    if(!is.matrix(l.gof$background)){
        l.gof$background=matrix(l.gof$background,ncol=length(l.gof$background),byrow=TRUE)
    }
    ## Compute the liklihood ratio per bin 
    trans.ll= apply(l.gof$transcribed,1,function(x) max(x)+log(sum(exp(x-max(x)))))
    untrans.ll= apply(l.gof$background,1,function(x) max(x)+log(sum(exp(x-max(x)))))
    lr=trans.ll-untrans.ll
    ## Check how many of the individual bin calls match the state of the most likely model
    match=numeric(length(lr))        
    match[which(lr>0)]=1
    trans=numeric(length(lr))        
    trans[indicator[[mod.map]]$trans]=1
    ## Get which invidiual bins matched the TU call bin states
    match.mod <- (match==trans)
    match.rle=rle(match.mod)
    return(list(match=sum(match.mod)/length(match.mod),lr=lr,match.rle=match.rle))
}

## Run TU select on a gene
geneTuSelect <-  function(ensg.id,txTable.gr,geneTileSums,trans.params,back.params,tile=args$tile){
    ## write(paste0("Selecting TU for ",ensg.id,"..."),stdout())
    zT=list(`+`="plus",`-`="minus")

    ## Pull out transcript models and bin sums relevant to the given ensg.id
    txMods=txTable.gr[txTable.gr$ensg.id==ensg.id,]
    zStr=zT[[strand(txMods)@values[1]]]
    txBinSums=geneTileSums[[zStr]][names(geneTileSums[[zStr]])==ensg.id,]
    ## Binarize binSums
    txBinSums$sum[txBinSums$sum>=1]=1
     ## Make line for untranscribed model
    utrans=GRanges(seqnames(txMods[1]),IRanges(min(start(txMods)),max(end(txMods)),names="Untranscribed"),
        strand=strand(txMods)[1],ensg.id=txMods$ensg.id[1],txtype=txMods$txtype[1],hngc.id=txMods$hngc.id[1])
    
    ## Compute the indicator function
    txInd=calc.indicator(txMods,txBinSums,tile)
    ## Reduce Models to identifiable set
    red=reduceModels(txMods,txInd)
    txMods=red$txMods
    txInd=red$txInd
    ## Augment models with potential untranscribed TU
    txInd[["untranscribed"]]=list(un.trans=seq(1,length(txBinSums)),trans=numeric(0))
    txMods=c(txMods,utrans)
    
    ## Compute posterior for each model 
    model.posterior=compute.model.posterior(trans.params,back.params,indicator=txInd,binned.data=txBinSums$sum)
    ## Select most likely model
    mod.map=as.integer(which.max(model.posterior))

    ## Compute goodness of fit statistics
    g.stat=gof.map(mod.map,txBinSums$sum,indicator=txInd)
    ## Compute longest non-matching run
    lnmr=max(c(g.stat$match.rle$lengths[g.stat$match.rle$values==FALSE],0))
    ## Convert models to output format, wierd formating changes for compatibility with previously
    ## written code.
    models=with(as.data.frame(txMods),data.table(TXCHROM=as.character(seqnames),TXSTART=start,TXEND=end,
        GENEID=ensg.id,TXNAME=names(txMods),TXSTRAND=as.character(strand),GENENAME=hngc.id,TXTYPE=txtype,
        post.prob=model.posterior[mod.map],gof=g.stat$match,lnmr=lnmr))

    ## Plot if appropriate -- REDOING PLOTTING WITH GVIZ AT LATER DATE
    ## if(geneid<=args$plot & plot.mod){
    ##    create.model.grob(create.models(models.aug,100),bins,model.posterior,indicator,g.stat,lnmr,as.character(models$GENENAME[1]),"",bin.size=100,buffer,antisense)
    ## }
    
    ## Return values, need to fix genename to hngc.id 
    return(models[mod.map])  
}

########################
## Functions for detecting runthrough
########################

## Splits TuList by strand and sorts it so correct detection of nearest upstream transcript can take place
splitSortTuList <-  function(tuList){
    ## split ...
    tuList.plus=tuList[TXSTRAND=="+"]
    tuList.minus=tuList[TXSTRAND=="-"]
    ## and sort
    tuList.plus=tuList.plus[with(tuList.plus,order(TXCHROM,TXSTART)),]
    tuList.minus=tuList.minus[with(tuList.minus,order(TXCHROM,-TXEND)),]
    ## and split again, this time on chromosomes
    tuList.plus.chr = splitDataTable(tuList.plus,"TXCHROM")
    tuList.minus.chr = splitDataTable(tuList.minus,"TXCHROM")
    return(list(p=tuList.plus.chr,m=tuList.minus.chr))
}

## Split data.table by a colummn
## http://stackoverflow.com/questions/14977997/split-data-table
splitDataTable <- function(dt,attr) {
  boundaries=c(0,which(head(dt[[attr]],-1)!=tail(dt[[attr]],-1)),nrow(dt))
  return(
    mapply(
      function(start,end) {dt[start:end,]},
      head(boundaries,-1)+1,
      tail(boundaries,-1),
      SIMPLIFY=F))
}

## Find nearest active transcript (NAT) to a gene
## Scan for NAT in nearest "search.range" previous transcripts
## If gene is returned as self-linked there is no NAT for that gene, otherwise
## The gene returned is the NAT and the LINK gene the gene it is the NAT of
findNAT <- function(subTuList,search.range=15){
    ## Reduce subTusList to only those that are transcribed
    subTuListTrans=subTuList[TXNAME!="Untranscribed"]
    ## Create a list to record 
    nat=subTuListTrans
    ## If chromosome has >=2 active transcripts
    if(nrow(subTuListTrans)>=2){
        ## If first element on list "add" to list with self-link
        for(i in 2:nrow(subTuListTrans)){
            nat[i]=nearest.active.transcript(i,subTuListTrans,search.range)
        }
        ## nat.dt=rbindlist(nat)
        ## Annotate NAT with the gene it is nearest to
        nat$LINK=subTuListTrans$GENEID
        ## Remove lines where there was no NAT
        ## nat.dt=nat.dt[GENEID!=LINK]
        return(nat)
    } else if(nrow(subTuListTrans)==1){
        ## When a chromosome only has one gene return it as self-linked
        nat$LINK=subTuListTrans$GENEID
        return(nat)
    }
}

## Get the nearest active transcript for
# gets the nearest active transcript on the same strand given a window range
# will also return transcripts that run through transcript of interest
# t.model.index is a line number from t.list which is the sorted df output
# Also only takes active transcripts
nearest.active.transcript <- function(index,t.sub.list,search.range){
    t.model=t.sub.list[index,]
    no.nat=t.model
    ## Generate search list
    t.search.list=t.sub.list[max(1,(index-search.range)):(index-1)]
    ## Now get the transcript with the end closest to the start of t.model
    if(t.model$TXSTRAND=="+"){
        ets=t.model$TXSTART-t.search.list$TXEND # compute end nat to start t.model tss 
        if(sum(ets>0)>0){
            nat.ind=which.min(ets[ets>0]) 
            nat=t.search.list[ets>0,][nat.ind,]
        } else {
            nat.ind=which.min(abs(ets)) 
            nat=t.search.list[nat.ind,]
        }
    } else {
        ets=t.search.list$TXSTART-t.model$TXEND # compute end nat to start t.model tss 
        if(sum(ets>0)>0){
            nat.ind=which.min(ets[ets>0]) 
            nat=t.search.list[ets>0,][nat.ind,]
        } else {
            nat.ind=which.min(abs(ets)) 
            nat=t.search.list[nat.ind,]
        }       
    }
    return(nat) 
}

## From the NAT and it's corresponding TU extract region between them
## Where dist(NAT,TUs)<=minDist, TU is added to auto-RT list
natToTu <- function(nat.p,nat.m,tuList,minDist=args$min_dist,maxDist=args$up_window){
    ## Match NATs to TUs
    natP=merge(tuList[,LINK:=GENEID],rbindlist(nat.p),by="LINK")
    natM=merge(tuList[,LINK:=GENEID],rbindlist(nat.m),by="LINK")

    ## Remove Genes where LINK==GENEID (aka. the first gene on each chromosome)
    pullnames=c("TXCHROM.x","TXSTART.x","TXEND.x","GENEID.x","TXSTRAND.x")
    auto.non.rt=rbind(natP[GENEID.y==LINK],natM[GENEID.y==LINK])[,pullnames,with=FALSE]
    setnames(auto.non.rt,c("chr","start","end","GENEID","strand"))
    natP=natP[GENEID.y!=LINK]
    natM=natM[GENEID.y!=LINK]
        
    ## Convert NATs and TUs to the respective RT regions
    ## TUs on the positive strand
    rt.p=with(natP,data.table(chr=TXCHROM.x,start=TXEND.y,end=TXSTART.x,GENEID=LINK,strand="+"))
    ## TUs on the negative strand
    rt.m=with(natM,data.table(chr=TXCHROM.x,start=TXEND.x,end=TXSTART.y,GENEID=LINK,strand="-"))

    ## Filter out genes where dist(NAT,TU)<=minDist and put them in auto-RT chunk of list
    auto.rt=rbind(rt.p[with(rt.p,end-start)<=minDist],rt.m[with(rt.m,end-start)<=minDist])
    rt.p=rt.p[!with(rt.p,end-start)<=minDist]
    rt.m=rt.m[!with(rt.m,end-start)<=minDist]

    ## Filter out genes where dist(NAT,TU)>=maxDist and put them in auto non-RT chunk of list
    auto.non.rt=rbind(auto.non.rt,rt.p[with(rt.p,end-start)>=maxDist],rt.m[with(rt.m,end-start)>=maxDist])
    rt.p=rt.p[!with(rt.p,end-start)>=maxDist]
    rt.m=rt.m[!with(rt.m,end-start)>=maxDist]
        
    return(list(p=rt.p,m=rt.m,auto=auto.rt,auto.non=auto.non.rt))
}

detect.runthrough <- function(rtCounts){
    ## Binarize data
    dat.m=rtCounts$sum
    dat.m[dat.m>=1]=1
    dat.m=data.frame(trans=dat.m)
    
    ##
    ## Specify two state model
    ##

    trans.mat=matrix(c(0.99,0.01,0.01,0.99),nrow=2)
    ## Emission prior
    emis.prior=c(0.6,0.1)
    state.prior=c(0.9,0.1)

    ## Specify which paramters will be fixed
    f=c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,FALSE,FALSE)
    ## Specify HMM model
    hmm = depmix(response=trans~1, data=dat.m,family = binomial(),nstates=nrow(trans.mat),instart=state.prior,trstart=trans.mat,respstart=emis.prior)
    ## Fit emission parameters
    invisible(hmm <- fit(hmm,fixed=f))
    ## Compute forward-backward table
    ## hmm.fb=forwardbackward(hmm)
    ## barplot(dat.m$trans)
    
    ##
    ## Specify one state model
    ##

    trans.mat.one=matrix(c(1),nrow=1)
    ## Estimate emission prior from start of rt-tu
    emis.prior.one=c(0.3)
    state.prior.one=c(1)
    
    ## Specify which paramters will be fixed
    f.one=c(TRUE,TRUE,FALSE)
    ## Specify HMM model
    hmm.one = depmix(response=trans~1, data=dat.m,family = binomial(),nstates=nrow(trans.mat.one),instart=state.prior.one,trstart=trans.mat.one,respstart=emis.prior.one)
    ## Fit emission parameters
    invisible(hmm.one <- fit(hmm.one,fixed=f.one))

    ## Compute the log-liklihood ratio test
    L01 <- as.vector(- 2 * (as.numeric(logLik(hmm.one)) - as.numeric(logLik(hmm))))
    p.val=pchisq(L01, df=1, lower.tail = FALSE)

    return(data.table(GENEID=as.character(names(rtCounts)[1]),rt.p=p.val))
}

## Given all tiling sums, split them into the correct number of chunks
splitRtTask <- function(rtTileSums,cores=args$parallel){
    temp=unlist(GRangesList(rtTileSums),use.names=FALSE)
    ## Allocate genes to cores
    core.id=data.table(gene=unique(names(temp)))
    core.id[,core.id:=rep(c(1:cores),length.out=length(gene))]
    setkey(x=core.id,"gene")

    ## Match genes with their cores in GRanges 
    cid=merge(data.table(id=1:length(temp),gene=names(temp)),core.id,by="gene")
    ## Maintain order of rows after merge
    cid=cid[order(id)]
    cid[,id:=NULL]
    temp$core.id=cid$core.id
    
    ## Split tile sums into list where each entry is a rt-tu
    temp <- split(temp,temp$core.id)
    return(temp)
}

## Do runthrough computations in parallel
parallel.rt <- function(reqData,cores=args$parallel){
    cl <- initParallel(cores)
    temp <- foreach(i=1:length(reqData),.combine="rbind") %dopar% {
        source("rtHMMLib.R")
        sub=reqData[[i]]
        rm(list="reqData")
        return(rbindlist(lapply(split(sub,names(sub)),function(x) detect.runthrough(x))))
    }
    stopCluster(cl)
    return(temp)
}
