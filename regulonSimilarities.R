




common_targets_n <- function(x, y = "SMARCA2", regulon)
length(intersect(regulon[regulon$tf == x,"target"], regulon[regulon$tf == y,"target"]))


regulon.w_2 <- addWeights(regulon = pruned.regulon_2,
      expMatrix = GeneExpressionMatrix[,selected_cells],
      exp_assay = "logcounts",
      peakMatrix = peakMatrix,
      peak_assay = "counts",
      exp_cutoff = NULL,
      clusters = GeneExpressionMatrix[,selected_cells]$Cellline,
      method = "wilcoxon")
regulon.w.motif_2 <- addMotifScore(regulon.w_2,
      archr_path=archR_project_path,
      species="human",
      genome="hg38")
regulon.w.moti_2f$weight[regulon.w.motif_2$motif==0,] <- 0
regulon.w.motif_2$weight[regulon.w.motif_2$motif==0,] <- 0
score.combine_2 <- calculateActivity(expMatrix = GeneExpressionMatrix[,selected_cells],
      regulon = regulon.w_2,
      mode = "weight",
      method = "weightedMean",
      exp_assay = "logcounts",
      normalize = FALSE)
cor_matr_2 <- cor(as.matrix(t(score.combine_2)))
head(sort(cor_matr_2["SMARCA2",], decreasing = TRUE), n =50)
cor_matr_2["SMARCA2","AR"]
hist(cor_matr_2["SMARCA2",])
cor_matr_2["SMARCA4","AR"]
head(sort(cor_matr_2["SMARCA4",], decreasing = TRUE), n =50)
length(intersect(regulon.w_2[regulon.w_2$tf == "SMARCA2","target"], regulon.w_2[regulon.w_2$tf == "SMC1A","target"]))
common_targets_n <- function(x, y = "SMARCA2", regulon)
length(intersect(regulon[regulon$tf == x,"target"], regulon[regulon$tf == y,"target"]))
unlist(lapply(sort(cor_matr_2["SMARCA2",], decreasing = TRUE)), function(x) common_targets_n(x=x, regulon = regulon.w_2))
unlist(lapply(sort(cor_matr_2["SMARCA2",], decreasing = TRUE), function(x) common_targets_n(x=x, regulon = regulon.w_2)))
res <- unlist(lapply(names(sort(cor_matr_2["SMARCA2",], decreasing = TRUE)), function(x) common_targets_n(x=x, regulon = regulon.w_2)))
res[1:10]
plot(res)
plot(res, pch=16)
plot(res, pch=16, cex = 0.5)
grep("FOX", rownames(score.combine), value = TRUE)
cor_matr_2["SMARCA4","FOXA1"]
cor_matr_2["SMARCA4","FOXA2"]
cor_matr_2["SMARCA4","FOXA3"]
cor_matr_2["SMARCA2","FOXA1"]
paste(names(head(sort(cor_matr_2["SMARCA4",], decreasing = TRUE), n=21)), as.character(unlist(lapply(names(head(sort(cor_matr_2["SMARCA4",], decreasing = TRUE), n=21)),
function(x) common_targets_n(x=x, regulon = regulon.w_2, y ="SMARCA4")))), sep = "=", collapse = ", ")
which(names(sort(cor_matr_2["SMARCA4",], decreasing = TRUE))=="AR")
which(names(sort(cor_matr_2["SMARCA4",], decreasing = TRUE))=="FOXA1")
