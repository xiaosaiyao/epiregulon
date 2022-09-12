#' A function to calculate transcriptional activity from a matrix of gene expression and regulon
#'
#' @param expMatrix a SummarizedExperiment object or matrix of gene expression with gene names matching those of the regulon or
#' @param regulon Matrix consisting of tf(regulator), target and a column indicating degree of association between TF and target
#                 such as "mor" or "corr".
#                 example regulon:
#                 tf      target  corr
#                 Esr1    Pgr     0.56
#' @param method String indicating the method used to calculate signature score. Available methods in gsva package include "gsva", "ssgsea", "zscore" and "plage"
#' @param mode String indicating the choice of TF-target association being either correlation "corr" or mode of regulation "mor".
#' Correlation can either be supplied or computed from expr using (corr_calculate = FALSE). Mor is supplied by the Dorothea regulon.
#' @param corr_calculate Logical indicating whether correlation needs to be calculated
#' @param save.regulon String indicating the path to save the new regulon with updated correlation
#' @param assay_name String indicating the name of the assay corresponding to gene expression
#' if expr is a SummarizedExperiment object
#' @param rowData_name String indicating the column name in the rowData that matches the gene names in the regulon
#' if expr is a SummarizedExperiment object
#' @param sample_n An integer indicating the number of samples to subsample if the original expression matrix is large
#'
#'
#' @return A matrix of inferred transcription factor (row) activities in samples (columns)
#' @export
#' @import utils GSVA SummarizedExperiment
#' @examples
#' # Load regulon from dorothea
#' library(dorothea)
#' data(dorothea_hs, package = "dorothea")
#' regulon <- dorothea_hs[which(dorothea_hs$confidence %in% c("A")),]
#'
#' # Create a mock expression matrix
#' genes <- unique(c(regulon$tf, regulon$target))
#' expr <- matrix(rnorm(length(genes)*20, mean=2), nrow = length(genes), ncol=20)
#' rownames(expr) <- genes
#' colnames(expr) <- paste("sample", 1:20)
#'
#' # Use mor embedded in Dorothea regulon
#' activity.mor <- calculateActivityBulk(expr, regulon, mode = "mor", method = "ssgsea")
#'
#' # Add correlation of Dorothea from user supplied gene matrix
#' activity.corr <- calculateActivityBulk(expr, regulon, mode = "corr",
#' method = "ssgsea", corr_calculate = TRUE)
#'
#' # Subsample
#' activity.corr <- calculateActivityBulk(expr, regulon, mode = "corr",
#' method = "ssgsea", corr_calculate = TRUE, sample_n = 10)
#'
#' # Use Dorothea regulon with pre-calculated calculation
#' \dontrun{
#' regulon <- readRDS(file = "/gstore/project/lineage/GRN/data/Dorothea_TCGA/regulon.Breast.TCGA.ABCDE.rds")
#' activity.corr <- calculateActivityBulk(expr, regulon, mode = "corr", method ="ssgsea")
#' }
#'
#' @description
#' This function generates the transcriptional activity from a matrix of bulk gene expression data and pre-existing regulons such as those
#' available from Dorothea. The direction of target genes (activating or repressive) can be defined either by the mode of regulation
#' included in the regulon or calculated based on the correlation between the expression of the TF and its target genes
#' @author Xiaosai Yao

calculateActivityBulk <- function(expMatrix,
                                  regulon,
                                  method = c("gsva", "ssgsea", "zscore", "plage"),
                                  mode = c("mor", "corr"),
                                  corr_calculate = NULL,
                                  save.regulon = NULL,
                                  sample_n = NULL,
                                  assay_name = "rpkm",
                                  rowData_name = "symbol") {

  # convert expr to matrix

    if (checkmate::test_class(expMatrix, classes = "SummarizedExperiment")){
      expr_assay <- assay(expMatrix, assay_name)
      if (!is.null(rowData_name)){
        rownames(expr_assay) <- rowData(expMatrix)[,rowData_name]
      }
      expr <- as.matrix(expr_assay)

    } else {
      expr <- as.matrix(expMatrix)
    }

    # remove tfs not found in expression matrix
    regulon <- regulon[(regulon$tf %in% rownames(expr)), ]

    # remove targets not found in expression matrix
    regulon <- regulon[(regulon$target %in% rownames(expr)), ]

    # split regulon by TF
    tf_indexes <- split(seq_len(nrow(regulon)), regulon$tf)
    unique_tfs <- names(tf_indexes)

    # add correlation
    if (mode == "corr") {
        regulon$corr = 0
        message("calculating correlation\n")
        pb <- txtProgressBar(min = 0,
                             max = length(unique_tfs),
                             style = 3)

        regulon_weight_list = vector("list", length(unique_tfs))
        names(regulon_weight_list) = unique_tfs

        counter = 0
        for (tf in unique_tfs) {

            if (is.null(sample_n) & (corr_calculate)) {
                tf_expr <- expr[tf, , drop = FALSE]
                target_expr_matrix <- expr[regulon$target[tf_indexes[[tf]]], , drop = FALSE]
            } else if (!is.null(sample_n) & (corr_calculate)) {
                samples <- sample(1:ncol(expr), sample_n, replace = FALSE)
                tf_expr <- expr[tf, samples , drop = FALSE]
                target_expr_matrix <- expr[regulon$target[tf_indexes[[tf]]], samples, drop = FALSE]
            }

            corr <- as.numeric(cor(t(tf_expr), t(target_expr_matrix), use = "everything"))
            regulon_weight_list[[tf]] = corr

            counter = counter + 1
            setTxtProgressBar(pb, counter)
            Sys.sleep(1 / 100)

        }
        regulon_weights <- unlist(regulon_weight_list)
        regulon$corr <- regulon_weights

        if (!is.null(save.regulon)){
            saveRDS(regulon, file = save.regulon)
        }

        regulon_pos <- regulon[which(regulon$corr > 0),]
        regulon_neg <- regulon[which(regulon$corr < 0),]
    }

    if (mode == "mor") {
        regulon_pos <- regulon[which(regulon$mor > 0),]
        regulon_neg <- regulon[which(regulon$mor < 0),]
    }



    regulon.ensg.pos.lst <- split(regulon_pos$target, regulon_pos$tf)
    regulon.ensg.neg.lst <- split(regulon_neg$target, regulon_neg$tf)

    #calculate positive and negative scores
    gsvaMatrix.pos <- gsva(
            expr,
            regulon.ensg.pos.lst,
            parallel.sz = 1,
            min.sz = 3,
            abs.ranking = F,
            method = method
        )
    gsvaMatrix.neg <- gsva(
            expr,
            regulon.ensg.neg.lst,
            parallel.sz = 1,
            min.sz = 3,
            abs.ranking = F,
            method = method
        )

    #merge positive and negative matrices
    all_genes <- S4Vectors::union(rownames(gsvaMatrix.pos), rownames(gsvaMatrix.neg))

    TFactivity.total <- matrix(0, ncol = ncol(gsvaMatrix.pos), nrow = length(all_genes))
    colnames(TFactivity.total) <- colnames(gsvaMatrix.pos)
    rownames(TFactivity.total) <- all_genes

    TFactivity.pos.total <- TFactivity.total
    TFactivity.neg.total <- TFactivity.total

    TFactivity.pos.total[match(rownames(gsvaMatrix.pos), rownames(TFactivity.total)),] <- gsvaMatrix.pos
    TFactivity.neg.total[match(rownames(gsvaMatrix.neg), rownames(TFactivity.total)),] <- gsvaMatrix.neg
    TFactivity.total <- TFactivity.pos.total - TFactivity.neg.total
}
