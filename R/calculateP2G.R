#' Establish peak to gene links based on correlations between ATAC-seq peaks and RNA-seq genes
#'
#' @param peakMatrix A SingleCellExperiment object containing counts of chromatin accessibility at each peak region or genomic bin from scATAC-seq.
#' `rowRanges` should contain genomic positions of the peaks in the form of `GRanges`.
#' @param expMatrix A SingleCellExperiment object containing gene expression counts from scRNA-seq. `rowRanges` should contain genomic positions of
#' the genes in the form of `GRanges`. `rowData` should contain a column of gene symbols with column name matching the `gene_symbol` argument.
#' @param reducedDim A matrix of dimension reduced values
#' @param cutoff_stat A names of a statistic used to determine significant links to assign peak to gene links.
#' Should be `Correlation`, `p_val` or `FDR`.
#' @param cutoff_sig A numeric scalar to specify the cutoff for the links between ATAC-seq peaks and RNA-seq genes .
#' Default is set to 0.5.
#' @param cellNum A numeric to specify the average number of cells per K-mean cluster. Alternatively, an object of the class `CellNumSol`
#' returned by `optimizeMetacellNumber` function. If set to `NULL`, its value is determined automatically, based on the number of cells.
#' @param maxDist An integer to specify the base pair extension from transcription start start for overlap with peak regions
#' @param exp_assay String indicating the name of the assay in expMatrix for gene expression
#' @param peak_assay String indicating the name of the assay in peakMatrix for chromatin accessibility
#' @param gene_symbol String indicating the column name in the rowData of expMatrix that corresponds to gene symbol
#' @param cor_method String indicating which correlation coefficient is to be computed. One of 'pearson' (default), 'kendall', or 'spearman'.
#' @param assignment_method String indicating the method used to assign target genes to regulatory elements. 'Correlation' is based on correlation between ATAC and RNA
#' above a correlation threshold set by cor_cutoff. 'Nearest' assigns the closest expressed gene to regulatory element meeting a correlation threshold set by cor_cutoff.
#' Set cor_cutoff to 0 if wishing to assign the closest expressed gene without any correlation cutoff
#' @param frac_RNA An integer to indicate the fraction of cells expressing a gene. It is used to filter the gene expression matrix for expressed genes
#' @param frac_ATAC An integer to indication the fraction of cells showing chromatin accessibility. It is used to filter the peak Matrix for open regions
#' @param nRandConns An integer specifying the number of false connections between regulatory elements and target genes which
#' will be used to calculate empirical p-values of correlation coefficients
#' @param BPPARAM A BiocParallelParam object specifying whether summation should be parallelized. Use BiocParallel::SerialParam() for
#' serial evaluation and use BiocParallel::MulticoreParam() for parallel evaluation
#' @param verbose A boolean indicating whether messages should be emitted during computation
#'
#' @return A DataFrame of Peak to Gene correlation
#' @importFrom SummarizedExperiment rowRanges rowData colData rowRanges<- rowData<-
#' @importFrom SingleCellExperiment altExp altExp<- reducedDim<-
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors Rle mcols mcols<- DataFrame
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @export
#'
#' @examples
#' # create a mock singleCellExperiment object for gene expression matrix
#' set.seed(1000)
#' gene_sce <- scuttle::mockSCE()
#' gene_sce <- scuttle::logNormCounts(gene_sce)
#' gene_gr <- GenomicRanges::GRanges(seqnames = Rle(c('chr1', 'chr2', 'chr3','chr4'), nrow(gene_sce)/4),
#'                    ranges = IRanges(start = seq(from = 1, length.out=nrow(gene_sce), by = 1000),
#'                    width = 100))
#' rownames(gene_sce) <- rownames(gene_sce)
#' gene_gr$name <- rownames(gene_sce)
#' rowRanges(gene_sce) <- gene_gr
#'
#' # create a mock singleCellExperiment object for peak matrix
#' peak_gr <- GenomicRanges::GRanges(seqnames = 'chr1',
#'                    ranges = IRanges(start = seq(from = 1, to = 10000, by = 1000), width = 100))
#' peak_counts <- matrix(sample(x = 0:4, size = ncol(gene_sce)*length(peak_gr), replace = TRUE),
#'                       nrow = length(peak_gr), ncol=ncol(gene_sce))
#' peak_sce <- SingleCellExperiment(list(counts = peak_counts), colData = colData(gene_sce))
#' rowRanges(peak_sce) <- peak_gr
#' rownames(peak_sce) <- paste0('peak',1:10)

#' # create a mock reducedDim matrix
#' reducedDim_mat <- matrix(runif(ncol(gene_sce)*50, min = 0, max = 1), nrow = ncol(gene_sce), 50)
#' p2g <- calculateP2G(peakMatrix = peak_sce, expMatrix = gene_sce, reducedDim = reducedDim_mat,
#'                     cellNum = 20)
#' @author Xiaosai Yao, Shang-yang Chen

calculateP2G <- function(peakMatrix = NULL,
                         expMatrix = NULL,
                         reducedDim = NULL,
                         cutoff_stat = c("p_val", "FDR", "Correlation"),
                         cutoff_sig = 0.05,
                         cellNum = NULL,
                         maxDist = 250000,
                         exp_assay = "logcounts",
                         peak_assay = "counts",
                         gene_symbol = "name",
                         cor_method = c("pearson", "spearman", "kendall"),
                         assignment_method = c("correlation","nearest"),
                         frac_RNA = 0,
                         frac_ATAC = 0,
                         nRandConns = 1e5,
                         BPPARAM = BiocParallel::SerialParam(progressbar = TRUE),
                         verbose = TRUE,
                         clusters=deprecated()
) {

    if (lifecycle::is_present(clusters)) {
        warning("Argument 'clusters' to calculateP2G was deprecated as of epiregulon version 2.0.0")
    }

    if(verbose){
        writeLines("Using epiregulon to compute peak to gene links...")
    }

    # check inputs
    cor_method <- match.arg(cor_method)
    assignment_method <- match.arg(assignment_method)
    cutoff_stat <- match.arg(cutoff_stat)
    .validate_input_sce(SCE=expMatrix, assay_name=exp_assay, row.ranges=TRUE)
    .validate_input_sce(SCE=peakMatrix, assay_name=peak_assay, row.ranges=TRUE)
    if (!identical(colnames(expMatrix), colnames(peakMatrix))){
        stop("Cell names in expMatrix and peakMatrix should be identical")
    }
    if(is.null(reducedDim)) stop("reducedDim argument is NULL.")

    if (!gene_symbol %in% colnames(rowData(expMatrix))) {
        stop("rowData of expMatrix does not contain ", gene_symbol)
    }
    if(class(cellNum)!="CellNumSol" && as.list(sys.call(sys.nframe()-1))[[1]]!="optimizeMetacellNumber"){
        message("Value of the paramater 'cellNum' has not been optimized.
                Consider running function 'optimizeMetacellNumber' and use output to set 'cellNum'")
    }
    if(class(cellNum)=="CellNumSol") {
        if (cellNum@args$cor_method != cor_method){
            warning(strwrap(sprintf("%s correlation method has been used for
                                     optimization of metacell number whereas
                                     `calculateP2G` is using %s method. It is
                                     recommended to use the same method in both
                                     functions.", cellNum@args$cor_method, cor_method)))
        }
        cellNum <- cellNum@solution^2
    }
    checkmate::check_double(cellNum,lower=1, upper=ncol(expMatrix))

    if(verbose){
        writeLines("Creating metacells...")
    }
    kNum = round(ncol(expMatrix)/cellNum)
    agg_data_list <- .create_metacells(expMatrix,
                                       exp_assay,
                                       peakMatrix,
                                       peak_assay,
                                       reducedDim,
                                       gene_symbol,
                                       frac_RNA,
                                       frac_ATAC,
                                       kNum=kNum)

    # find overlaps between REs and resized TGs
    if(verbose){
        writeLines("Looking for regulatory elements near target genes...")
    }
    o <- .find_nearby_REs(geneStart=agg_data_list[["geneStart"]],
                          peakSet=agg_data_list[["peakSet"]],
                          expMatrix_agg=agg_data_list[["geneExpr"]],
                          peakMatrix_agg=agg_data_list[["peakCounts"]],
                          old.idxRNA=agg_data_list[["old.idxRNA"]],
                          old.idxATAC=agg_data_list[["old.idxATAC"]],
                          maxDist=maxDist,
                          gene_symbol=gene_symbol,
                          assignment_method=assignment_method)

    # Calculate correlation
    if(verbose){
        writeLines("Computing correlations...")
    }

    o$Correlation <- NA

    idx_pairs <- mapply(function(x,y) list(c(x,y)), as.integer(o$RNA), as.integer(o$ATAC))
    split_points <- seq(1,length(idx_pairs), by = 2e4)
    o$Correlation <- unlist(BiocParallel::bplapply(X = split_points,
                                        FUN = .RE_TG_correlation,
                                        idx_pairs=idx_pairs,
                                        exprMatrix=agg_data_list[["geneExpr"]],
                                        peakMatrix=agg_data_list[["peakCounts"]],
                                        cor_method=cor_method,
                                        BPPARAM = BPPARAM))

    o <- .addFDR(df=o,
                geneStart = agg_data_list[["geneStart"]],
                peakSet=agg_data_list[["peakSet"]],
                geneExpr = agg_data_list[["geneExpr"]],
                peakCounts = agg_data_list[["peakCounts"]],
                n_random_conns = nRandConns,
                cor_method = cor_method,
                BPPARAM = BPPARAM)

    p2g_merged <- o[, c("old.idxATAC", "chr", "start", "end", "old.idxRNA", "Gene","Correlation", "p_val", "FDR", "distance")]
    colnames(p2g_merged) <- c("idxATAC", "chr", "start", "end", "idxRNA", "target","Correlation", "p_val", "FDR", "distance")
    if(cutoff_stat=="Correlation"){
        relation_fun <- get(">")
    }
    else{
        relation_fun <- get("<")
    }
    p2g_merged <- p2g_merged[relation_fun(p2g_merged[,cutoff_stat], cutoff_sig), , drop = FALSE]

    p2g_merged <- p2g_merged[order(p2g_merged$idxATAC, p2g_merged$idxRNA), , drop = FALSE]
    return(p2g_merged)
}


#' @importFrom GenomicRanges resize mcols rowRanges
#' @importFrom scrapper aggregateAcrossCells clusterKmeans
.create_metacells <- function(expMatrix, exp_assay, peakMatrix, peak_assay, reducedDim,
                              gene_symbol, frac_RNA, frac_ATAC, kNum){

    kclusters <- clusterKmeans(t(as.matrix(reducedDim)),k = kNum)$clusters
    kclusters <- as.character(kclusters)
    geneStart <- resize(rowRanges(expMatrix), width=1)
    mcols(geneStart)[,gene_symbol] <- rowData(expMatrix)[,gene_symbol]
    data_to_aggregate <- as(assay(expMatrix, exp_assay), "CsparseMatrix")
    # aggregate by k-means clusters
    res <- aggregateAcrossCells(data_to_aggregate, factors = list(kclusters))
    expMatrix <- t(t(res$sums)/res$counts)

    peakSet = rowRanges(peakMatrix)
    data_to_aggregate <- as(assay(peakMatrix, peak_assay), "CsparseMatrix")
    res <- aggregateAcrossCells(data_to_aggregate, factors = list(kclusters))
    peakMatrix <- t(t(res$sums)/res$counts)
    # keep track of the original ATAC and expression indices
    old.idxRNA <- seq_len(nrow(expMatrix))
    old.idxATAC <- seq_len(nrow(peakMatrix))

    # filter gene expression matrix based on fraction of cells expressing gene
    cells_express_rna <- rowSums(expMatrix > 0)
    frac_expressed_rna <- cells_express_rna/ncol(expMatrix)
    expMatrix <- expMatrix[frac_expressed_rna > frac_RNA, ]
    old.idxRNA <- old.idxRNA[frac_expressed_rna > frac_RNA]
    geneStart <- geneStart[frac_expressed_rna > frac_RNA]

    # filter peak matrix based on fraction of cells showing chromatin accessibility
    cells_express_atac <- rowSums(peakMatrix > 0)
    frac_expressed_atac <- cells_express_atac/ncol(peakMatrix)
    peakMatrix <- peakMatrix[frac_expressed_atac > frac_ATAC, ]
    old.idxATAC <- old.idxATAC[frac_expressed_atac > frac_ATAC]
    peakSet <- peakSet[frac_expressed_atac > frac_ATAC]

    # return gene expression and peak matrix
    return(list(geneExpr = expMatrix,
                peakCounts = peakMatrix,
                geneStart=geneStart,
                peakSet = peakSet,
                old.idxRNA = old.idxRNA,
                old.idxATAC = old.idxATAC))
}

#' @importFrom GenomicRanges findOverlaps start end distance distanceToNearest mcols

.find_nearby_REs <- function(geneStart, peakSet, expMatrix_agg, peakMatrix_agg,
                             old.idxRNA, old.idxATAC,
                             maxDist, gene_symbol, assignment_method){
    if (assignment_method == "correlation"){
        o <- DataFrame(findOverlaps(resize(geneStart, maxDist, "center"),
                                    peakSet, ignore.strand = TRUE))
    } else if (assignment_method == "nearest") {

        # assign every regulatory element to its nearest gene
        nearest_gene <- DataFrame(distanceToNearest(peakSet, geneStart))

        # flip order of query and subject to be consistent with correlation mode
        o <- DataFrame(queryHits=nearest_gene$subjectHits,
                       subjectHits=nearest_gene$queryHits,
                       distance=nearest_gene$distance)

        # filter by distance
        o <- o[which(o$distance < maxDist), ]

    }

    #Get Distance from Fixed point A B
    o$distance <- distance(geneStart[o[, 1]], peakSet[o[, 2]])
    colnames(o) <- c("RNA", "ATAC", "distance")

    # add old idxRNA and idxATAC
    o$old.idxRNA <- old.idxRNA[o[, 1]]
    o$old.idxATAC <- old.idxATAC[o[, 2]]

    #add metadata to o
    o$Gene <- mcols(geneStart)[o[, 1], gene_symbol]
    o$chr <- as.character(seqnames(peakSet)[o[, 2]])
    o$start <- start(peakSet)[o[, 2]]
    o$end <- end(peakSet)[o[, 2]]
    return(o)
}
#' @importFrom GenomicRanges seqnames

.addFDR <- function(df,
                   geneStart,
                   peakSet,
                   geneExpr,
                   peakCounts,
                   n_random_conns,
                   cor_method,
                   BPPARAM){
    # take a sample from the marginal distribution of peaks in RE-TG connections
    random_peak_idx <- sample(df[,"ATAC"], n_random_conns, replace=TRUE)
    seq_peaks <- unique(seqnames(peakSet[random_peak_idx]))
    aligned_random_peaks <- c()
    aligned_random_genes <- c()
    for(seq_peak in seq_peaks){
        seq_peak_idx <- which(as.logical(seqnames(peakSet[random_peak_idx])==seq_peak))
        # find genes in other chromosomes
        remote_gene_idx <- which(as.logical(seqnames(geneStart)!=seq_peak))
        aligned_random_genes <- c(aligned_random_genes, sample(remote_gene_idx, length(seq_peak_idx), replace=TRUE))
        aligned_random_peaks <- c(aligned_random_peaks, random_peak_idx[seq_peak_idx])
    }
    # tie matching genes and peaks into pairs
    idx_pairs <- mapply(function(x,y) list(c(x,y)), aligned_random_genes, aligned_random_peaks)
    # determine chunk limits for parallelization
    split_points <- seq(1,length(idx_pairs), by = 2e4)
    null_correlations <- unlist(BiocParallel::bplapply(X = split_points,
                                                      FUN = .RE_TG_correlation,
                                                      idx_pairs=idx_pairs,
                                                      exprMatrix=geneExpr,
                                                      peakMatrix=peakCounts,
                                                      cor_method=cor_method,
                                                      BPPARAM = BPPARAM))

    rand_corr_distr_pos <- ecdf(null_correlations[null_correlations>=0])
    rand_corr_distr_neg <- ecdf(null_correlations[null_correlations<=0])
    correlations <- df$Correlation
    df$p_val <- rep(NA, nrow(df))
    df$p_val[correlations==0] <- 1
    df$p_val[correlations<0] <- rand_corr_distr_neg(correlations[correlations<0])
    df$p_val[correlations>0] <- (1-rand_corr_distr_pos(correlations[correlations>0]))
    df$FDR <- NA
    df$FDR[df$Correlation==0] <- 1
    df$FDR[df$Correlation<0] <- p.adjust(df$p_val[df$Correlation<0], method="BH")
    df$FDR[df$Correlation>0] <- p.adjust(df$p_val[df$Correlation>0], method="BH")
    return(df)
}

#' Determine the optimal number of metacells to be used by `calculateP2G` function
#'
#' This function attempts to find the optimal value for the `cellNum` parameter, which
#' is used by `calculateP2G` to define the number of metacells. The value of this
#' parameter is critical: too many clusters may lead to insufficient signal integration
#' to overcome the effect of data sparsity, whereas too few clusters may result in
#' excessive averaging and loss of important biological variability.
#'
#' @param peakMatrix A SingleCellExperiment object containing counts of chromatin accessibility at each peak region or genomic bin from scATAC-seq.
#' `rowRanges` should contain genomic positions of the peaks in the form of `GRanges`.
#' @param expMatrix A SingleCellExperiment object containing gene expression counts from scRNA-seq. `rowRanges` should contain genomic positions of
#' the genes in the form of `GRanges`. `rowData` should contain a column of gene symbols with column name matching the `gene_symbol` argument.
#' @param reducedDim A matrix of dimension reduced values
#' @param exp_assay String indicating the name of the assay in expMatrix for gene expression
#' @param peak_assay String indicating the name of the assay in peakMatrix for chromatin accessibility
#' @param subsample_prop A numeric indicating the fraction of features from `expMatrix` and `peakMatrix` used to optimize `kNum`
#' @param n_iter An integer indicating the number of iterations before in which the value of `kNum` parameter is optimized
#' @param cellNumMin A numeric used to optimize value of `cellNum` parameter. Corresponds to the lower bound for the
#' average number of cells per K-mean cluster in the first iteration of the optimization algorithm. If `cellNum` is not `NULL`
#' this parameter is ignored.
#' @param cellNumMax A numeric used to optimize value of `cellNum` parameter. Corresponds to the upper bound for the
#' average number of cells per K-mean cluster in the first iteration of the optimization algorithm. If `cellNum` is not `NULL`
#' this parameter is ignored.
#' @param n_evaluation_points An integer defining how many metacells numbers are tested in the first iteration to find
#' the optimal one. Must not be less than 3. If `n_inter` > 1 new evaluation points (metacell numbers) are added in the proximity of the current solution.
#' @param ... Other arguments passed to `calculateP2G` function
#'
#' @return An object of the class `CellNumSol` to be passed to `calculateP2G` as `cellNum` paramater.
#' @export

optimizeMetacellNumber <- function(peakMatrix,
                                      expMatrix,
                                      reducedDim,
                                      exp_assay,
                                      peak_assay,
                                      subsample_prop=0.1,
                                      n_iter=1,
                                      cellNumMin=NULL,
                                      cellNumMax=NULL,
                                      n_evaluation_points=4,
                                      ...){
    # check inputs
    .validate_input_sce(SCE=expMatrix, assay_name=exp_assay, row.ranges=TRUE)
    .validate_input_sce(SCE=peakMatrix, assay_name=peak_assay, row.ranges=TRUE)
    if (!identical(colnames(expMatrix), colnames(peakMatrix))){
        stop("Cell names in expMatrix and peakMatrix should be identical")
    }
    if(is.null(reducedDim)) stop("reducedDim argument is NULL.")
    checkmate::check_double(subsample_prop,lower=0, upper=1)
    checkmate::check_integer(n_evaluation_points, lower=3)

    n_cells <- ncol(peakMatrix)
    if(is.null(cellNumMin)){
        cells_per_cluster_min <- min(20, round(n_cells/10))
    }
    else{
        cells_per_cluster_min <- cellNumMin
    }

    if(is.null(cellNumMax)){
        cells_per_cluster_max <- min(2000, round(n_cells/10))
    }
    else{
        cells_per_cluster_max <- min(cellNumMax, n_cells/3)
    }
    cells_per_cluster_min <- min(cells_per_cluster_min, cells_per_cluster_max)

    assay(expMatrix, exp_assay) <- as(assay(expMatrix, exp_assay), "CsparseMatrix")
    assay(peakMatrix, peak_assay) <- as(assay(peakMatrix, peak_assay), "CsparseMatrix")

    if(subsample_prop<1){
        selected_peak_idx <- sort(sample(nrow(peakMatrix), round(subsample_prop*nrow(peakMatrix))))
        peakMatrix <- peakMatrix[selected_peak_idx,]
        selected_gene_idx <- sort(sample(nrow(expMatrix), round(subsample_prop*nrow(expMatrix))))
        expMatrix <- expMatrix[selected_gene_idx,]
    }

    evaluation_points <- seq(sqrt(cells_per_cluster_min), sqrt(cells_per_cluster_max),
                             length.out=n_evaluation_points)
    # drop evaluation points that are duplicates after mapping to cluster numbers
    kNum <- round(n_cells/evaluation_points^2)
    evaluation_points <- evaluation_points[!duplicated(kNum)]
    if(length(evaluation_points)<3){
        stop("To few evaluation points to optimize kNum paramater. Consider using more cells or changing cellNumMin or cellNumMax parameters.")
    }
    areas <- c()
    for (i in seq_along(evaluation_points)){
        p2g <- calculateP2G(
            peakMatrix = peakMatrix,
            expMatrix = expMatrix,
            reducedDim = reducedDim,
            exp_assay = exp_assay,
            peak_assay = peak_assay,
            cellNum = evaluation_points[i]^2,
            verbose = FALSE,
            cutoff_sig = 2,
            ...
        )
        p_val_pos_reg <- p2g$p_val[p2g$Correlation>=0] # should (all) zeros be included?
        # calculate area under p-value cumulative distribution curve
        areas <- c(areas, mean(p_val_pos_reg))
    }
    regr_data = data.frame(areas=areas, sqrt_cellNum=evaluation_points)
    lin_model <- lm(areas~poly(sqrt_cellNum,2,raw=TRUE), data=regr_data)
    sol <- optim(17, function(x) predict(lin_model,
                                          newdata=data.frame(sqrt_cellNum=x)),
                 method="Brent",
                 lower=min(regr_data$sqrt_cellNum),
                 upper=max(regr_data$sqrt_cellNum))$par
    last_iteration = 1L
    if(n_iter>1){
        for(n in 2:n_iter){
            # select three additional evaluation points from the interval
            # containing current solution and two adjacent ones
            sol_position <- findInterval(sol, evaluation_points)
            # extend domain to account for solution being in the first or last interval
            extended_domain <- c(1,evaluation_points, sqrt(round(n_cells/3)))
            # get the points delimiting intervals that will be split
            # (interval with solution and adjacent ones)
            interval_limits <- sol_position + 0:3
            evaluation_points_new <- extended_domain[interval_limits]
            # calculate mid-points of the selected intervals
            evaluation_points_new <- evaluation_points_new[-length(evaluation_points_new)] + diff(evaluation_points_new)/2
            kNum_new <- round(n_cells/(evaluation_points_new^2))
            # find kNum values that have already been checked
            already_used_filter <- kNum_new %in% round(n_cells/evaluation_points^2)
            kNum_new <- kNum_new[!already_used_filter]
            if(length(kNum_new)==0) break
            # adjust evaluation points to kNum
            evaluation_points_new <- evaluation_points_new[!already_used_filter]
            # duplicates might be generated as a result of rounding
            evaluation_points_new[!duplicated(kNum_new)]
            areas_new <- c()
            # calculate AUC (mean p-value) for each evaluation point
            for (i in seq_along(evaluation_points_new)){
                p2g <- calculateP2G(
                    peakMatrix = peakMatrix,
                    expMatrix = expMatrix,
                    reducedDim = reducedDim,
                    exp_assay = exp_assay,
                    peak_assay = peak_assay,
                    cellNum = evaluation_points_new[i]^2,
                    verbose = FALSE,
                    cutoff_sig = 2,
                    ...
                )
                p_val_pos_reg <- p2g$p_val[p2g$Correlation>=0]
                # calculate area under p-value cumulative distribution curve
                areas_new <- c(areas_new, mean(p_val_pos_reg))
            }
            evaluation_points <- c(evaluation_points, evaluation_points_new)
            areas <- c(areas, areas_new)
            areas <- areas[order(evaluation_points)]
            evaluation_points <- sort(evaluation_points)
            regr_data = data.frame(areas=areas, sqrt_cellNum=evaluation_points)
            lin_model <- lm(areas~poly(sqrt_cellNum, 2, raw=TRUE), data=regr_data)
            sol <- optim(300, function(x) predict(lin_model,
                                                newdata=data.frame(sqrt_cellNum=x)),
                         method="Brent",
                         lower=min(regr_data$sqrt_cellNum),
                         upper=max(regr_data$sqrt_cellNum))$par
            last_iteration <- last_iteration+1L
        }
    }
    # get values of some arguments used in the function call
    # as that might help in troubleshooting
    arg_names <- c("n_iter", "n_evaluation_points", "cellNumMin",
              "cellNumMax", "subsample_prop")
    default_args <- arg_names[!arg_names %in% names(as.list(match.call()))]
    args <- as.list(match.call())[setdiff(arg_names, default_args)]
    # add default arguments from the function definition
    args <- c(args, formals(sys.function())[default_args])
    args <- args[arg_names]
    p2g_args = formals(calculateP2G)[c("nRandConns", "cor_method", "maxDist", "frac_RNA", "frac_ATAC")]
    user_specified_args <- intersect(names(p2g_args), names(list(...)))
    # replace deafults with the user specified arguments passed to calculateP2G
    p2g_args[user_specified_args] <- list(...)[user_specified_args]
    estimation_issue <- FALSE
    if(lin_model$coefficients[3] <= 0){
        warning("Coefficient of quadratic term in linear regression is not potitive.")
        estimation_issue <- TRUE
    }
    if(any(abs(sol-range(evaluation_points))<1e-4)){
        warning("Solution at the boundary of examined range.")
        estimation_issue <- TRUE
    }
    if(summary(lin_model)$r.squared<0.7) { # check if points are not too far from the curve
        warning("Coefficient of determination of the regression model is lower thant 0.7")
        estimation_issue <- TRUE
    }
    if(estimation_issue){
        # TO DO: reference to the on-line documentation
        message(paste(c(strwrap("An issue detected during estimation optimal number of metacells.
                        Consider at least one of the following actions:"),
                        "1. Change of the `cellNumMin` and `cellNumMax` paramaters",
                        "2. Increasing the number of evaluation points (`n_evaluation_points` argument)",
                        "3. Increasing the number of iterations (`n_iter` argument)",
                        strwrap("4. Increasing the number of false connections used to compute
                        p-value null distribution (`nRandConns` argument)"))))
        message("Solution not found using quadratic regression. Using cluster size with the lowest mean p-value.")
        sol <- evaluation_points[which.min(areas)]
    }

    new("CellNumSol", solution=sol,
        evaluation_points=evaluation_points,
        AUC = areas,
        regr_coefficients = lin_model$coefficients,
        n_cells = n_cells,
        last_iteration=last_iteration,
        args=c(args, p2g_args)
        )
}

.RE_TG_correlation <- function(ind, idx_pairs, exprMatrix, peakMatrix, cor_method){
    # select pairs to be included in this batch
    idx_pairs <- idx_pairs[ind:min(ind+2e4-1,length(idx_pairs))]
    vapply(idx_pairs,
           function(idx_pair) stats::cor(exprMatrix[idx_pair[1],],
                                                                 peakMatrix[idx_pair[2],],
                                                                 method=cor_method),
           numeric(1))
}

setMethod("plot", signature=c(x="CellNumSol"), function(x){
    x_val <- x@evaluation_points
    y_val <- x@AUC
    df <- data.frame(x=x_val, y=y_val)
    df2 <- data.frame(x=seq(min(x_val), max(x_val), length.out=1e3))
    df2$y <- cbind(1, df2$x, df2$x^2) %*% x@regr_coefficients
    plot(y~x, data=df2, type="l", xlab="Square root of the number of cells per cluster",
             ylab="Area under curve", ylim=range(c(df2$y,y_val)))
    sol=x@solution
    points(x_val,y_val,pch=16)
    lines(c(sol, sol), range(c(df2$y,y_val)), lt=2, col="red")
})



