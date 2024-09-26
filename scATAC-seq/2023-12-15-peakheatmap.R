#BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")

library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
peak <- readPeakFile("/home/xi/Desktop/aggregated_scATAC_summits.bed")
covplot(peak, chr = c("chr1", "chr2"))

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
# 定义TSS上下游的距离
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(peak, windows=promoter)
tagHeatmap(tagMatrix,xlim = c(-3000,3000),color = "red")

peakHeatmap(peak, TxDb=txdb, 
            upstream=3000, downstream=3000)
