edge_weights <- seq(0,1,length.out=20)
edge_weights[5:12] <- rev(edge_weights[5:12])
edge_weights[9] <- 0
edge_weights[7] <- NA

regulon = data.frame(tf = rep(LETTERS[1:10], each = 2), idxATAC = rep(1:10, 2), target  = rep(LETTERS[10:14], 4), corr = edge_weights)
regulon$target[10:13] <- LETTERS[15]
# duplicated tf-re and tf-tg pairs for aggregation
regulon[regulon$tf == "B", "idxATAC"] <- 3
regulon[regulon$tf == "C", "idxATAC"] <- 5
regulon[regulon$tf == "I", "target"] <- "K"
regulon[regulon$tf == "J", "target"] <- "M"

test_graph <- buildGraph(regulon, mode = "tripartite", weights = "corr")

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

test_that("buildGraph function correctly connects vertices in in the tripartite mode", {
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
  re_list[[tf]] <- split(re_list[[tf]],re_list[[tf]]$idxATAC) |>
    lapply(function(x) data.frame(idxATAC = x$idxATAC[1], corr = x$corr[which.max(abs(x$corr))])) |>
    {function(x) do.call(rbind, x)}()
  re_list[[tf]] <- re_list[[tf]][base::order(re_list[[tf]]$corr),,drop = FALSE]
  re_list[[tf]] <- re_list[[tf]][base::order(re_list[[tf]]$idxATAC),,drop = FALSE]
  rownames(re_list[[tf]]) <- NULL
}

for(idxATAC in sort(base::unique(regulon$idxATAC))){
  target_list[[idxATAC]]<- regulon[regulon$idxATAC == idxATAC, c("target", "corr"),drop = FALSE]
  target_list[[idxATAC]][is.na(target_list[[idxATAC]][,"corr"]),"corr"] <- 0
  # aggregate edges
  target_list[[idxATAC]] <- split(target_list[[idxATAC]],target_list[[idxATAC]]$target) |>
    lapply(function(x) data.frame(target = x$target[1], corr = x$corr[which.max(abs(x$corr))])) |>
    {function(x) do.call(rbind, x)}()
  target_list[[idxATAC]] <- target_list[[idxATAC]][base::order(target_list[[idxATAC]]$corr),,FALSE]
  target_list[[idxATAC]] <- target_list[[idxATAC]][base::order(target_list[[idxATAC]]$target),,FALSE]
  rownames(target_list[[idxATAC]]) <- NULL
}

test_graph <- buildGraph(regulon, mode = "tripartite", weights = "corr")

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

test_that("buildGraph function correctly assigns weights in the tripartite mode", {
  expect_identical(re_list_graph, re_list)
  expect_identical(target_list_graph, target_list)
})

regulon_na_rem <- regulon
regulon_na_rem[is.na(regulon_na_rem$corr),"corr"] <- 0

# create second data frame defining network structure
regulon_2 <- regulon

# change part of the weights to have non-zero edges after graphs subtraction
regulon_2[regulon_2$target %in% c("K", "O", "J", "L"), "corr"] <-
  seq(0,1, length.out = sum(regulon_2$target %in% c("K", "O", "J", "L")))

# for some node part of the edge should be 0 after graph subtraction
regulon_2[which(regulon_2$target=="L")[1], "corr"] <- regulon[which(regulon$target=="L")[1], "corr"]
regulon_2[which(regulon_2$target=="L")[1:2], "corr"] <- regulon[which(regulon$target=="O")[1:2], "corr"]

aggregate_table <- function(tab, vertices_cols){
  pairs <- paste0(tab[[vertices_cols[1]]], "_", tab[[vertices_cols[2]]])
  split(tab, pairs) |> lapply(function(x) data.frame(c1 = x[[vertices_cols[1]]][1],
                                                     c2 = x[[vertices_cols[2]]][1],
                                                     corr = x$corr[which.max(abs(x$corr))])) |>
    {function(x) do.call(rbind, x)}() -> tab
  colnames(tab) <- c(vertices_cols, "corr")
  tab
}

# aggregate separately tf-re and re-tg pairs and calculate weight difrrerences
tf_re_table_1 <- regulon_na_rem[,c("tf", "idxATAC", "corr")] |> aggregate_table(c("tf","idxATAC"))
tf_re_table_2 <- regulon_2[,c("tf", "idxATAC", "corr")] |> aggregate_table(c("tf","idxATAC"))
tf_re_diff_table <- tf_re_table_1[,c("tf", "idxATAC")]
tf_re_diff_table$corr <- abs(tf_re_table_1$corr - tf_re_table_2$corr)
tf_re_diff_table <- tf_re_diff_table[tf_re_diff_table$corr!=0,]

re_tg_table_1 <- regulon_na_rem[,c("idxATAC", "target", "corr")] |> aggregate_table(c("idxATAC", "target"))
re_tg_table_2 <- regulon_2[,c("idxATAC", "target", "corr")] |> aggregate_table(c("idxATAC", "target"))
re_tg_diff_table <- re_tg_table_1[,c("idxATAC", "target")]
re_tg_diff_table$corr <- abs(re_tg_table_1$corr - re_tg_table_2$corr)
re_tg_diff_table <- re_tg_diff_table[re_tg_diff_table$corr!=0,]
diff_table <- merge(tf_re_diff_table[,1:2], re_tg_diff_table[1:2])

diff_regulon <- data.frame(tf=regulon$tf, idxATAC = regulon$idxATAC, target = regulon$target, corr=abs(regulon_na_rem$corr - regulon_2$corr))
diff_regulon <- diff_regulon[diff_regulon$corr!=0,]

test_graph_1 <- buildGraph(regulon, mode = "tripartite", weights = "corr")
test_graph_2 <- buildGraph(regulon_2, mode = "tripartite", weights = "corr")
test_diff_graph <- buildDiffGraph(test_graph_1, test_graph_2)

# find re and tg associated with each transcription factor
re_list <- split(diff_table$idxATAC, diff_table$tf) |> lapply(unique) |> lapply(sort) |> lapply(as.character)
target_list <- split(diff_table$target, diff_table$tf) |> lapply(unique) |> lapply(sort)

# the same based on the graph
re_list_graph <- list()
target_list_graph <- list()
tfs <- igraph::V(test_diff_graph)[igraph::V(test_diff_graph)$type == "transcription factor"]$name
for(tf in sort(base::unique(tfs))){
  nearby_nodes <- igraph::ego(test_diff_graph, nodes = igraph::V(test_diff_graph)[igraph::V(test_diff_graph)$type == "transcription factor" & igraph::V(test_diff_graph)$name == tf], order = 2, mode = "out")
  re_list_graph[[tf]] <- sort(as.character(base::unique(nearby_nodes[[1]][nearby_nodes[[1]]$type == "regulatory element"]$name)))
  target_list_graph[[tf]] <- sort(as.character(base::unique(nearby_nodes[[1]][nearby_nodes[[1]]$type == "target gene"]$name)))
}

test_that("buildDiffGraph function correctly connects vertices in in the tripartite mode", {
  expect_identical(re_list_graph, re_list)
  expect_identical(target_list_graph, target_list)
})

# test builDiffGraph with respect to weight assignment
# find regulatory elements associated with respective tfs, and target genes associated with regulatory elements
re_list <- split(tf_re_diff_table, tf_re_diff_table$tf) |>
  lapply(function(x) x[order(x$corr),c("idxATAC", "corr")]) |>
  lapply(function(x) {x <- x[order(x$idxATAC),]; rownames(x) <- NULL; x})

target_list <- split(re_tg_diff_table, re_tg_diff_table$idxATAC) |>
  lapply(function(x) x[order(x$corr),c("target", "corr")]) |>
  lapply(function(x) {x <- x[order(x$target),]; rownames(x) <- NULL; x})

test_graph <- buildDiffGraph(test_graph_1, test_graph_2)

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

test_that("buildDiffGraph function correctly assigns weights in the tripartite mode", {
  expect_identical(re_list_graph, re_list)
  expect_identical(target_list_graph, target_list)
})
