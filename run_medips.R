# this part is run on a big server with lots of disk/memory

library("BSgenome")
library("MEDIPS")
library(BSgenome.Hsapiens.UCSC.hg19)
BSgenome="BSgenome.Hsapiens.UCSC.hg19"
uniq=1e-3
extend=300
shift=0
ws=100
# modify this
#chr.select="chr22"

# list all files
bams <- list.files(".",pattern="bam$")

# read em all in
xx <- sapply(bams,function(x) { 
  MEDIPS.createSet(file = x,
  BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq,
  window_size = ws)
})

names(xx) <- gsub(".bam","",names(xx))

# need to read in the sample groups based on "medical help to conceive"
s <- read.table("samplegroups.txt")

# separate two groups
mset1 <- xx[which(names(xx) %in% s[which(s$V2=="no"),1] )]
mset2 <- xx[which(names(xx) %in% s[which(s$V2=="yes"),1] )]

## CBMC first - select "C" samples first
mset1 <- mset1[grep("C",names(mset1))]
mset2 <- mset2[grep("C",names(mset2))]

# prepare the coupling set - not sure what that is but what the hey
CS <- MEDIPS.couplingVector(pattern = "CG", refObj = xx[[1]])

# run the contrast
mr.edgeR <- MEDIPS.meth(
  MSet1 = mset1, 
  MSet2 = mset2,
  CSet = CS, 
  p.adj = "fdr",
  diff.method = "edgeR",
  MeDIP = TRUE,
  CNV = FALSE,
  minRowSum = 10)

mr.edgeR.c <- mr.edgeR[order(mr.edgeR$edgeR.p.value),]

write.table(mr.edgeR.c, file="mr.edgeR.c.tsv",sep="\t",quote=FALSE)

# WBC
# separate two groups
mset1 <- xx[which(names(xx) %in% s[which(s$V2=="no"),1] )]
mset2 <- xx[which(names(xx) %in% s[which(s$V2=="yes"),1] )]

## CBMC first - select "C" samples first
mset1 <- mset1[grep("W",names(mset1))]
mset2 <- mset2[grep("W",names(mset2))]

# prepare the coupling set - not sure what that is but what the hey
CS <- MEDIPS.couplingVector(pattern = "CG", refObj = xx[[1]])

# run the contrast
mr.edgeR <- MEDIPS.meth(
  MSet1 = mset1,
  MSet2 = mset2,
  CSet = CS,
  p.adj = "fdr",
  diff.method = "edgeR",
  MeDIP = TRUE,
  CNV = FALSE,
  minRowSum = 10)

mr.edgeR.w <- mr.edgeR[order(mr.edgeR$edgeR.p.value),]

write.table(mr.edgeR.w, file="mr.edgeR.w.tsv",sep="\t",quote=FALSE)

save.image("castillo.Rdata")
