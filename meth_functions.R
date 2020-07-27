suppressPackageStartupMessages({
  library("plyr")
  library("R.utils")
  library("missMethyl")
  library("limma")
  library("DMRcate")
  library("topconfects")
  library("minfi")
  library("IlluminaHumanMethylation450kmanifest")
  library("RColorBrewer")
  library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
  library("eulerr")
  library("plyr")
  library("gplots")
  library("reshape2")
  library("beeswarm")
  library("RCircos")
})


# scree plot shows the amount of variation in a dataset that is accounted
# for by the first N principal components
myscree <- function(mx,n=10,main="") {
  pc <- princomp(mx)$sdev
  pcp <- pc/sum(pc)*100
  pcp <- pcp[1:10]
  barplot(pcp,cex.names = 1,las=2,ylim=c(0,60),
      ylab="percent (%) variance explained", main=main)
  text((0.5:length(pcp)*1.2),pcp,label=signif(pcp,3),pos=3,cex=0.8)
}


# Here is a function to make a volcano plot
make_volcano <- function(dm,name,mx) {
    sig <- subset(dm,adj.P.Val<0.05)
    N_SIG=nrow(sig)
    N_UP=nrow(subset(sig,logFC>0))
    N_DN=nrow(subset(sig,logFC<0))
    HEADER=paste(N_SIG,"@5%FDR,", N_UP, "up", N_DN, "dn")
    plot(dm$logFC,-log10(dm$P.Val),cex=0.5,pch=19,col="darkgray",
        main=name, xlab="log FC", ylab="-log10 pval")
    mtext(HEADER)
    grid()
    points(sig$logFC,-log10(sig$P.Val),cex=0.5,pch=19,col="red")
}

# Here is a function to make heatmaps 
make_heatmap <- function(dm,name,mx,n) {
  topgenes <-  rownames(head(dm[order(dm$P.Value),],n))
  ss <- mx[which(rownames(mx) %in% topgenes),]
  my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 25)
  heatmap.2(ss,scale="row",margin=c(10, 10),cexRow=0.4,trace="none",cexCol=0.4,
      col=my_palette, main=name)
}

# make beeswarm charts
# dm = a limma differential meth object
# name = character name of the limma dm object
# mx = matrix of normalised data
# groups = a vector of factors corresponding to the cols in mx
# n = the number of top significant genes to plot (default = 15) 
make_beeswarms <- function(dm,name,mx,groups,n=15) {
    par(mar=c(3,3,1,1))
    NCOLS=5
    NROWS=floor(n/NCOLS)
    if (n %% NCOLS > 0) { NROWS <- NROWS + 1 }
    par(mfrow=c(NROWS, NCOLS))
    topgenes <-  rownames(head(dm[order(dm$P.Value),],n))
    ss <- mx[which(rownames(mx) %in% topgenes),]
    n <- 1:n
    g1name=levels(groups)[1]
    g2name=levels(groups)[2]
    g1dat <- ss[n,which(groups == g1name)]
    g2dat <- ss[n,which(groups == g2name)]
    g1l <-lapply(split(g1dat, row.names(g1dat)), unlist)
    g2l <-lapply(split(g2dat, row.names(g2dat)), unlist)

    for (i in n) {
      mydat <- list(g1l[[i]],g2l[[i]])
        beeswarm(mydat,ylim=c(0,1),cex=0.2, pch=19,
        las=2, cex.lab=0.6, main=names( g1l )[i] , 
        ylab="",labels = c(g1name,g2name))
      grid()
    }
}

# make beeswarm charts for best confects
# dm = a limma differential meth object
# name = character name of the limma dm object
# mx = matrix of normalised data
# groups = a vector of factors corresponding to the cols in mx
# n = the number of top significant genes to plot (default = 15) 
make_beeswarms_confects <- function(confects,name,mx,groups,n=15) {
    par(mar=c(3,3,1,1))
    NCOLS=5
    NROWS=floor(n/NCOLS)
    if (n %% NCOLS > 0) { NROWS <- NROWS + 1 }
    par(mfrow=c(NROWS, NCOLS))
    topgenes <-  head(confects$table,n)$name
    ss <- mx[which(rownames(mx) %in% topgenes),]
    n <- 1:n
    g1name=levels(groups)[1]
    g2name=levels(groups)[2]
    g1dat <- ss[n,which(groups == g1name)]
    g2dat <- ss[n,which(groups == g2name)]
    g1l <-lapply(split(g1dat, row.names(g1dat)), unlist)
    g2l <-lapply(split(g2dat, row.names(g2dat)), unlist)

    for (i in n) {
      mydat <- list(g1l[[i]],g2l[[i]])
        beeswarm(mydat,ylim=c(0,1),cex=0.2, pch=19,
        las=2, cex.lab=0.6, main=names( g1l )[i] , 
        ylab="",labels = c(g1name,g2name))
      grid()
    }
}

# this is a wrapper which creates three charts
# We will be adding more
make_dm_plots <- function(dm,name,mx,groups=groups,confects=confects) {
    make_volcano(dm,name,mx)
    make_beeswarms(dm ,name , mx , groups , n= 15)
    make_heatmap(dm , name , mx ,n = 50)
    make_beeswarms_confects(confects, name, mx, groups, n=15)
}  

# this function performs DMRcate for peak calling
run_dmrcate <- function(mx=mxs,design=design) {
    fit.reduced <- lmFit(mxs,design)
    fit.reduced <- eBayes(fit.reduced)

    dmg <- makeGenomicRatioSetFromMatrix(mxs, rownames = NULL, pData = NULL ,
        array = "IlluminaHumanMethylation450k" ,
        mergeManifest = FALSE, what = "M")
    
    myannotation <- cpg.annotate("array", dmg, arraytype = "450K",
        analysis.type="differential", design=design, coef=3)

    dmrcoutput <- dmrcate(myannotation, lambda=1000, C=3)

    dmr <- extractRanges(dmrcoutput, genome = "hg19")

    return(dmr)
}


# this is a function which will perform differential methylation analysis
# if you provide it with the right inputs
dm_analysis <- function(samplesheet,sex,groups,mx,name,myann,beta) {

    design <- model.matrix(~ sex + groups)
    mxs <- mx[,which( colnames(mx) %in% samplesheet$Basename )]
    fit.reduced <- lmFit(mxs,design)
    fit.reduced <- eBayes(fit.reduced)
    summary(decideTests(fit.reduced))
    dm <- topTable(fit.reduced,coef=3, number = Inf)
    dma <- merge(myann,dm,by=0)
    dma <- dma[order(dma$P.Value),]
    dm_up <- rownames(subset(dm,adj.P.Val<0.05 & logFC>0))
    dm_dn <- rownames(subset(dm,adj.P.Val<0.05 & logFC<0))
    length(dm_up)
    length(dm_dn)
    confects <- limma_confects(fit.reduced, coef=3, fdr=0.05)
    dmr <- run_dmrcate(mx=mxs,design=design) 
    head(dmr)
    make_dm_plots(dm = dm ,name=name , mx=beta, groups= groups, confects=confects)
    dat <- list("dma"=dma, "dm_up"=dm_up, "dm_dn"=dm_dn, "confects"=confects, "dmr"= dmr)
    return(dat)
}

# Spearman ranks
myranks <- function(x) {
  xx<-x[[1]]
  xx$score <- sign(xx$logFC)/log10(xx$adj.P.Val)
  y <- xx[,"score",drop=FALSE]
  y$rn <- xx$Row.names
  return(y)
}

# enrichment analysis functions
run_mitch_rank <-function(dma){
  dmap <- dma[grep("Promoter_Associated",dma$Regulatory_Feature_Group),]
  dmap[which(dmap$UCSC_RefGene_Name==""),2] <- "NA"
  dmap$genename <- sapply(strsplit(dmap$UCSC_RefGene_Name,";"),"[[",1)
  dmap2 <- dmap[,c("genename","t")]
  #rank <- aggregate(. ~ genename,dmap2,function(x) { 
  #  x[which.max(abs(x))]
  #  } )
  rank <- aggregate(. ~ genename,dmap2,mean)
  rownames(rank) <- rank$genename
  rank$genename=NULL
  return(rank)
}

run_mitch_1d <- function(dma,name) {
  library("mitch")
  rank <- run_mitch_rank(dma)
  capture.output(
    res <- mitch_calc(x = rank,genesets = genesets, priority = "significance",resrows=20)
    , file = "/dev/null", append = FALSE,
    type = c("output", "message"), split = FALSE)
  head(res$enrichment_result,20)
  capture.output(
    mitch_plots(res,outfile=paste(name,".pdf"))
    , file = "/dev/null", append = FALSE,
    type = c("output", "message"), split = FALSE)
  return(res$enrichment_result)
}