edge_weights <- seq(0,1,length.out=20)
edge_weights[5:12] <- rev(edge_weights[5:12])
edge_weights[9] <- 0
edge_weights[7] <- NA

regulon = data.frame(tf = rep(LETTERS[1:10], each = 2), idxATAC = rep(1:10, 2), target  = rep(LETTERS[10:14], 4), corr = edge_weights)
regulon$target[10:13] <- LETTERS[15]
test_graph <- buildGraph(regulon)

re_list <- list()
target_list <- list()

# determine which regulatory elements and target genes can be reached from each transcription factor
# note that graph will disassemble original tf-re-tg triplets into tf-re and re-tg pairs.
for(tf in sort(base::unique(regulon$tf))){
  re_list[[tf]] <- sort(base::unique(as.character(regulon[regulon$tf == tf, "idxATAC"])))
  target_list[[tf]] <- sort(base::unique(regulon[regulon$idxATAC %in% re_list[[tf]], "target"]))
}

re_list_graph <- list()
target_list_graph <- list()
tfs <- igraph::V(test_graph)[igraph::V(test_graph)$type == "transcription factor"]$name
for(tf in sort(base::unique(tfs))){
  nearby_nodes <- igraph::ego(test_graph, nodes = igraph::V(test_graph)[igraph::V(test_graph)$type == "transcription factor" & igraph::V(test_graph)$name == tf], order = 2, mode = "out")
  re_list_graph[[tf]] <- sort(as.character(base::unique(nearby_nodes[[1]][nearby_nodes[[1]]$type == "regulatory element"]$name)))
  target_list_graph[[tf]] <- sort(as.character(base::unique(nearby_nodes[[1]][nearby_nodes[[1]]$type == "target gene"]$name)))
}

test_that("buildGraph function correctly connect vertices", {
  expect_identical(re_list_graph, re_list)
  expect_identical(target_list_graph, target_list)
})

# find regulatory elements associated with respective tfs and target genes associated with regulatory elements
re_list <- list()
target_list <- list()
regulon$idxATAC <- as.character(regulon$idxATAC)
for(tf in sort(base::unique(regulon$tf))){
  re_list[[tf]] <- regulon[regulon$tf == tf, c("idxATAC", "corr"),drop = FALSE]
  re_list[[tf]][is.na(re_list[[tf]][["corr"]]),"corr"] <- 0
  # aggregate tf-re pairs
  re_list[[tf]] <- split(re_list[[tf]], re_list[[tf]]$idxATAC) |> lapply(function(x) data.frame(idxATAC = x$idxATAC[1], corr = mean(x$corr)))
  re_list[[tf]] <- do.call(rbind, re_list[[tf]])
  re_list[[tf]] <- re_list[[tf]][base::order(re_list[[tf]]$corr),,drop = FALSE]
  re_list[[tf]] <- re_list[[tf]][base::order(re_list[[tf]]$idxATAC),,drop = FALSE]
  rownames(re_list[[tf]]) <- NULL
}

for(idxATAC in sort(base::unique(regulon$idxATAC))){
  target_list[[idxATAC]]<- regulon[regulon$idxATAC == idxATAC, c("target", "corr"),drop = FALSE]
  target_list[[idxATAC]][is.na(target_list[[idxATAC]][,"corr"]),"corr"] <- 0
  # aggregate re-tg pairs
  target_list[[idxATAC]] <- split(target_list[[idxATAC]], target_list[[idxATAC]]$target) |> lapply(function(x) data.frame(target = x$target[1], corr = mean(x$corr)))
  target_list[[idxATAC]] <- do.call(rbind, target_list[[idxATAC]])
  target_list[[idxATAC]] <- target_list[[idxATAC]][base::order(target_list[[idxATAC]]$corr),,FALSE]
  target_list[[idxATAC]] <- target_list[[idxATAC]][base::order(target_list[[idxATAC]]$target),,FALSE]
  rownames(target_list[[idxATAC]]) <- NULL
}

test_graph <- buildGraph(regulon)

tfs <- igraph::V(test_graph)[igraph::V(test_graph)$type == "transcription factor"]$name
re_list_graph <- list()
edges <- igraph::as_edgelist(test_graph)
for(tf in sort(base::unique(tfs))){
  # find neighbouring nodes of tf of interest
  nearby_nodes <- igraph::ego(test_graph, nodes = igraph::V(test_graph)[igraph::V(test_graph)$type == "transcription factor" & igraph::V(test_graph)$name == tf],
                              order = 1, mode = "out", mindist = 1)[[1]]
  from_ind <- which(edges[,1]==tf)
  df <- data.frame()
  # get weights of the edges  the neighbouring nodes
  for(i in seq_along(nearby_nodes)){
    to_ind <- which(edges[,2] == nearby_nodes[i]$name)
    neighbour_weights = E(test_graph)[intersect(from_ind, to_ind)]$weight
    df <- rbind(df, data.frame("idxATAC" = nearby_nodes[i]$name, "corr" = neighbour_weights))
  }
  df <- df[base::order(df$corr),]
  re_list_graph[[tf]] <- df[base::order(df$idxATAC),]
  rownames(re_list_graph[[tf]]) <- NULL
}

reg_elems <- igraph::V(test_graph)[igraph::V(test_graph)$type == "regulatory element"]$name
target_list_graph <- list()
for(re in sort(base::unique(reg_elems))){
  # find neighbouring downstream nodes of re of interest
  nearby_nodes <- igraph::ego(test_graph, nodes = igraph::V(test_graph)[igraph::V(test_graph)$type == "regulatory element" & igraph::V(test_graph)$name == re],
                              order = 1, mode = "out", mindist = 1)[[1]]
  from_ind <- which(edges[,1]==re)
  df <- data.frame()
  # get weights of the edges  the neighbouring nodes
  for(i in seq_along(nearby_nodes)){
    to_ind <- which(edges[,2] == nearby_nodes[i]$name)
    neighbour_weights = E(test_graph)[intersect(from_ind, to_ind)]$weight
    df <- rbind(df, data.frame("target" = nearby_nodes[i]$name, "corr" = neighbour_weights))
  }
  df <- df[base::order(df$corr),]
  target_list_graph[[re]] <- df[base::order(df$target),]
  rownames(target_list_graph[[re]]) <- NULL
}

test_that("buildGraph function correctly assign weights in tripartite mode", {
  expect_identical(re_list_graph, re_list)
  expect_identical(target_list_graph, target_list)
})
