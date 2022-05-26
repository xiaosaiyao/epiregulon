#' A function to calculate weights for the regulons by computing co-association between TF and target gene expression
#'
#' @param regulon A data frame consisting of tf (regulator) and target in the column names. Additional columns indicating degree
#' of association between tf and target such as "mor" or "corr" are optional.
#' @param sce A SingleCellExperiment object containing gene expression information
#' @param cluster_factor String specifying the field in the colData of the SingleCellExperiment object to be averaged as pseudobulk (such as cluster)
#' @param block_factor String specifying the field in the colData of the SingleCellExperiment object to be used as blocking factor (such as batch)
#' @param exprs_values String specifying the name of the assay to be retrieved from the SingleCellExperiment object
#' @param corr Logical scalar indicating whether to calculate weights based on correlation
#' @param MI Logical scalar indicating whether to calculate weights based on mutual information
#' @param min_targets Integer specifying the minimum number of targets for each tf in the regulon with 10 targets as the default
#' @param BPPARAM A BiocParallelParam object specifying whether summation should be parallelized. Use BiocParallel::SerialParam() for
#' serial evaluation and use BiocParallel::MulticoreParam() for parallel evaluation
#' @param alt.exp A matrix that is used in place of gene expression to correlate with target gene expression. See details.
#' @param alt.exp.merge A logical to indicate whether to consider both TF expression and alt.exp matrix. See details.
#'
#' @return A data frame with columns of corr and/or MI added to the regulon. TFs not found in the expression matrix and regulons not
#' meeting the minimal number of targets were filtered out.
#' @importFrom stats cor
#' @details
#' The default mode is to compute weights by correlating the pseudobulk target gene expression vs the pseudobulk TF gene expression.
#' However, often times, an inhibitor of TF does not alter the gene expression of the TF. In rare cases, cells may even compensate
#' by increasing the expression of the TF. In this case, the activity of the TF, if computed by gene expression correlation, may show a
#' spuriously increase. As an alternative to gene expression, we may use accessibility associated with TF, such as those computed by
#' chromVar. When alt.exp.merge is true, we take the product of the gene expression and the values in the alt.exp matrix.
#'
#'
#' @export
#'
#' @examples
#' # create a mock singleCellExperiment object for gene expression matrix
#' example_sce <- scuttle::mockSCE()
#' example_sce <- scuttle::logNormCounts(example_sce)
#' example_sce$cluster <- sample(LETTERS[1:5], ncol(example_sce), replace = TRUE)
#'
#' # create a mock regulon
#' regulon <- data.frame(tf = c(rep("Gene_0001",5), rep("Gene_0002",10)),
#'                       target = c(paste0("Gene_000",2:6), paste0("Gene_00",11:20)))
#'
#' # add weights to regulon
#' regulon.w <- addWeights(regulon, example_sce, cluster_factor="cluster", exprs_values = "logcounts",
#'                         min_targets = 5)
#'
#' # Alternatively, add a matrix of chromVar values in place of TF expression
#' \dontrun{
#' #create chromVar values from archR
#' library(ArchR)
#' library(parallel)
#' proj <- addBgdPeaks(proj)
#' proj <- addDeviationsMatrix(ArchRProj = proj,
#' peakAnnotation = "motif", force = TRUE,logFile = "addDeviation")
#'
#' #retrieve chromVar matrix
#' #chromVar
#' motifMatrix= getMatrixFromProject(
#' ArchRProj = proj,
#' useMatrix = "motifMatrix",
#' useSeqnames = NULL,
#' verbose = TRUE,
#' binarize = FALSE,
#' threads = getArchRThreads(),
#' logFile = createLogFile("getMatrixFromProject")
#' )
#'
#' # calculate weights using alt.exp
#' regulon.w.2 <- addWeights(regulon, example_sce, cluster_factor = "cluster",
#' exprs_values = "logcounts", min_targets = 5, alt.exp = assay(motifMatrix, "z"),
#' alt.exp.merge = TRUE)
#'}

addWeights = function(regulon,
                      sce,
                      cluster_factor,
                      block_factor = NULL,
                      exprs_values = "logcounts",
                      corr = TRUE,
                      MI = FALSE,
                      min_targets = 10,
                      BPPARAM = BiocParallel::SerialParam(),
                      alt.exp = NULL,
                      alt.exp.merge = FALSE){


  # define groupings
  groupings = DataFrame(cluster = colData(sce)[cluster_factor])
  if (!is.null(block_factor)) {
    groupings$block = colData(sce)[block_factor]
  }

  # compute average expression across clusters and batches
  writeLines("calculating average expression across clusters...")

  averages.se = scater::sumCountsAcrossCells(
    sce,
    exprs_values = exprs_values,
    ids = groupings,
    average = T,
    BPPARAM = BPPARAM
  )


  if (!is.null(alt.exp)) {
    alt.avg.se = scater::sumCountsAcrossCells(
      alt.exp,
      ids = groupings,
      average = T,
      BPPARAM = BPPARAM
    )
    alt.avg = assays(alt.avg.se)$average
  }
  # average expression across pseudobulk clusters
  expr = assays(averages.se)$average

  # remove genes whose expressions are NA for all pseudobulks
  expr = expr[!rowSums(is.na(expr)) == ncol(expr), ]

  # order regulon
  regulon = regulon[order(regulon$tf, regulon$target),]

  # remove tfs not found in expression matrix or alt.exp matrix
  if (is.null(alt.exp)){
    regulon = regulon[(regulon$tf %in% rownames(expr)),]
  } else {
    regulon = regulon[(regulon$tf %in% rownames(alt.exp)),]
  }

  if (alt.exp.merge) {
    regulon = regulon[(regulon$tf %in% intersect(rownames(expr), rownames(alt.exp))),]
  }

  # remove tfs with less than min_targets
  regulon = regulon[regulon$tf %in% names(which(table(regulon$tf) >= min_targets)),]

  #remove targets not found in expression matrix
  regulon = subset(regulon, (target %in% rownames(expr)))

  tf_indexes = split(seq_len(nrow(regulon)), regulon$tf)
  unique_tfs = names(tf_indexes)


  # compute correlation
  if (corr) {
    writeLines("computing correlation of the regulon...")
    pb = txtProgressBar(min = 0,
                        max = length(unique_tfs),
                        style = 3)
    counter = 0

    regulon_weight_list = vector("list", length(unique_tfs))
    names(regulon_weight_list) = unique_tfs


    for (tf in unique_tfs) {
      if (is.null(alt.exp)) {
        tf_expr = expr[tf, ,drop = FALSE ]
      } else {
        tf_expr = alt.avg[tf, , drop=FALSE]
      }

      if (alt.exp.merge){
        tf_expr = expr[tf, ,drop = FALSE ]
        tf_alt = alt.avg[tf, , drop=FALSE]
        tf_expr = tf_expr * tf_alt
      }
      target_expr_matrix = expr[regulon$target[tf_indexes[[tf]]], ,drop = FALSE ]
      weights = as.numeric(cor(t(tf_expr), t(target_expr_matrix), use = "everything"))
      regulon_weight_list[[tf]] = weights
      Sys.sleep(1 / 100)
      counter = counter + 1
      setTxtProgressBar(pb, counter)

    }
    regulon_weights = unlist(regulon_weight_list)
    regulon$weight = regulon_weights

  }
  # compute mutual information
  if (MI) {
    writeLines("computing mutual information of the regulon...")

    n_pseudobulk =  length(unique(colData(sce)[,cluster_factor]))

    if (n_pseudobulk < 5) {
      stop("Too few clusters for mutual information calculation. Need at least 5 clusters")

    } else{

      regulon_MI = c()
      counter = 0

      pb = txtProgressBar(min = 0,
                          max = length(unique_tfs),
                          style = 3)

      regulon_MI_list = vector("list", length(unique_tfs))
      names(regulon_MI_list) = unique_tfs

      for (tf in unique_tfs) {

        if (is.null(alt.exp)) {
          tf_expr = expr[tf, ,drop = FALSE ]
        } else {
          tf_expr = alt.avg[tf, , drop=FALSE]
        }

        target_expr_matrix = expr[regulon$target[tf_indexes[[tf]]], ]


        if (length(unique(expr[tf,])) <  5) {
          MI = rep(NA,nrow(target_expr_matrix))

        }else{
          MI = numeric(nrow(target_expr_matrix))

          for (i in seq_len(nrow(target_expr_matrix))) {
            target = rownames(target_expr_matrix)[i]
            if (length(unique(expr[target,])) <  5) {
              MI[i] = NA
            } else{
              y2d = entropy::discretize2d(expr[tf,],
                                          expr[target,],
                                          numBins1 =  n_pseudobulk,
                                          numBins2 =  n_pseudobulk)
              MI[i] = entropy::mi.empirical(y2d)

            }
          }


        }

        regulon_MI_list[[tf]] = MI


        Sys.sleep(1 / 100)
        counter = counter + 1
        setTxtProgressBar(pb, counter)

      }
      regulon_MI = unlist(regulon_MI_list)
      regulon$MI = regulon_MI
    }

  }


  return(regulon)

}
