aggregateMatrix <- function(regulon, mode, FUN) {

    groupings <- paste(regulon$tf, regulon$target, sep = "#")
    groupings <- factor(groupings, levels = unique(groupings))

    weights <- regulon[, mode]
    agg.weights <- rowsum(weights, groupings, reorder = FALSE)

    if (FUN == "mean") {
        num <- table(groupings)
        agg.weights <- agg.weights/as.integer(num[rownames(agg.weights)])
    }


    rownames.split <- do.call(rbind, strsplit(rownames(agg.weights),
        "#"))
    aggregated <- S4Vectors::DataFrame(tf = rownames.split[, 1],
        target = rownames.split[, 2], weight = I(agg.weights))
    colnames(aggregated)[3] <- mode
    aggregated
}

renameCluster <- function(clusters) {
    if (!is.null(clusters)) {
        clusters[clusters == "all"] <- "clusters_all"
    }
    clusters
}

initiateMatCluster <- function(clusters, nrow, value = NA) {
    unique_clusters <- sort(unique(clusters))
    cluster_mat <- matrix(value, nrow = nrow, ncol = length(unique_clusters) +
        1)
    colnames(cluster_mat) <- c("all", unique_clusters)
    cluster_mat
}





