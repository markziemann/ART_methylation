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
chr.select="chr22"

# list all files
bams <- list.files(".",pattern="bam$")

# read em all in
xx <- sapply(bams,function(x) { 
  MEDIPS.createSet(file = x,
  BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq,
  window_size = ws)
})

# separate two groups
mset1 <- xx[grepl("C",names(xx))]
mset2 <- xx[grepl("W",names(xx))]

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

mr.edgeR <- mr.edgeR[order(mr.edgeR$edgeR.p.value),]

save.image("castillo.Rdata")
