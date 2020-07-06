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

bams <- list.files(".",pattern="bam")

xx <- sapply(bams,function(x) { 
  MEDIPS.createSet(file = x,
  BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq,
  window_size = ws, chr.select = chr.select)
})

mr.edgeR = MEDIPS.meth(
  MSet1 = , 
  MSet2 = ,
  CSet = CS, 
  p.adj = "bonferroni",
  diff.method = "edgeR",
  MeDIP = TRUE,
  CNV = F,
  minRowSum = 10)
