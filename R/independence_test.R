test_independence <- function(regulon, PeakMatrix, GeneExpressionMatrix, cluster = "C5"){
    cluster_cells <- which(colData(GeneExpressionMatrix)$Clusters == cluster)
    GeneExpressionMatrix <- assay(GeneExpressionMatrix)[,cluster_cells]
    PeakMatrix <- assay(PeakMatrix)[,cluster_cells]
    range_list <- list()
    i<-0
    while(i*1000 < nrow(regulon)){
        range_list[[i+1]] <- c(1000*(i)+1, min(1000*(i+1), nrow(regulon)))
        i <- i +1
    }
    counts <- BiocParallel::bplapply(X = range_list,
                                     FUN = test_independence_bp,
                                     regulon,
                                     PeakMatrix,
                                     GeneExpressionMatrix,
                                     cluster_cells,
                                     cluster,
                                    BPPARAM = BiocParallel::MulticoreParam(workers=8))


    counts <- do.call(rbind, counts)
    counts <- counts[order(counts$i),]
    counts <- counts[,colnames(counts)!="i"]
    regulon <- cbind(regulon, counts)
    regulon$joint_prob_pred <- regulon$tf_re * regulon$tg
    regulon
}



test_independence_bp <- function(range_vector, regulon, PeakMatrix, GeneExpressionMatrix, cluster_cells, cluster = "C5"){
    counts <- data.frame()
    for (i in range_vector[1]:range_vector[2]){
        tf_active <- GeneExpressionMatrix[regulon[i, "tf"],]>1
        tg_active <- GeneExpressionMatrix[regulon[i, "idxRNA"], ]>1
        re_available <- PeakMatrix[regulon[i, "idxATAC"], ] > 0
        counts <- rbind(counts, data.frame(i = i, tf_re = sum((tf_active+re_available)==2)/length(cluster_cells),
                   tg = sum(tg_active)/length(cluster_cells),
                   triple = sum((tf_active+tg_active+re_available)==3)/length(cluster_cells)))
    }
    counts
}


reprogram_seq_data_joint_prob <- test_independence(reprogram_seq_data, PeakMatrix, GeneExpressionMatrix)



add_p_val <- function(regulon, n_cells){
    p_vals<-c()
    for (i in seq_len(nrow(regulon))){
        if (regulon$triple[i] == 0){
            p_vals[length(p_vals)+1] <- 1
            next
        }
        p_vals[length(p_vals)+1] <- 1 - pbinom(n_cells*regulon[i,"triple"]-1, n_cells, regulon[i,"joint_prob_pred"])
    }
    regulon$p_val <- p_vals
    regulon$p_val_adj <- regulon$p_val*n_cells
    regulon$p_val_adj[regulon$p_val_adj>1] <- 1
    regulon
}

n <- length(which(colData(GeneExpressionMatrix)$Clusters=="C5"))

reprogram_seq_data_joint_prob <- add_p_val(reprogram_seq_data_joint_prob, n)

colnames(reprogram_seq_data_joint_prob)[colnames(reprogram_seq_data_joint_prob) == "p_val"] <- paste0("p_val_", "C5")
colnames(reprogram_seq_data_joint_prob)[colnames(reprogram_seq_data_joint_prob) == "p_val_adj"] <- paste0("p_val_adj_", "C5")
colnames(reprogram_seq_data_joint_prob)[colnames(reprogram_seq_data_joint_prob) == "triple"] <- paste0("triple_", "C5")
colnames(reprogram_seq_data_joint_prob)[colnames(reprogram_seq_data_joint_prob) == "joint_prob_pred"] <- paste0("joint_prob_pred_", "C5")
colnames(reprogram_seq_data_joint_prob)[colnames(reprogram_seq_data_joint_prob) == "tg"] <- paste0("tg_", "C5")
colnames(reprogram_seq_data_joint_prob)[colnames(reprogram_seq_data_joint_prob) == "tf_re"] <- paste0("tf_re_", "C5")

add_new_cluster <- function(cluster, old_regulon){
    regulon <- test_independence(old_regulon, PeakMatrix, GeneExpressionMatrix, cluster = cluster)
    n = length(which(colData(GeneExpressionMatrix)$Clusters==cluster))
    regulon <- add_p_val(regulon, n)
    colnames(regulon)[colnames(regulon) == "p_val"] <- paste0("p_val_", cluster)
    colnames(regulon)[colnames(regulon) == "p_val_adj"] <- paste0("p_val_adj_", cluster)
    colnames(regulon)[colnames(regulon) == "triple"] <- paste0("triple_", cluster)
    colnames(regulon)[colnames(regulon) == "joint_prob_pred"] <- paste0("joint_prob_pred_", cluster)
    colnames(regulon)[colnames(regulon) == "tg"] <- paste0("tg_", cluster)
    colnames(regulon)[colnames(regulon) == "tf_re"] <- paste0("tf_re_", cluster)
    regulon
}

reprogram_seq_data_joint_prob <- add_new_cluster("C1", reprogram_seq_data_joint_prob)
reprogram_seq_data_joint_prob <- add_new_cluster("C3", reprogram_seq_data_joint_prob)


