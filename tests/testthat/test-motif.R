regulon <- S4Vectors::DataFrame(tf = c('AR','AR','AR','ESR1','ESR1','NKX2-1'),
idxATAC = 1:6)
peaks <- GenomicRanges::GRanges(seqnames = c('chr12','chr19','chr19','chr11','chr6','chr1'),
ranges = IRanges::IRanges(start = c(124914563,50850845, 50850844, 101034172, 151616327, 1000),
end = c(124914662,50850929, 50850929, 101034277, 151616394,2000)))
regulon <- addMotifScore(regulon, peaks=peaks)

test_that("addMotifScore works", {
  expect_equal(regulon$motif, c(0, 0, 0, 1, 1, 0))
})
