regulon <- data.frame(tf = rep(LETTERS[1:4], times = c(5,5,6,3)), target = NA)
regulon$target[regulon$tf=="A"] <- LETTERS[5:9]
regulon$target[regulon$tf=="B"][c(1,3,5)] <- regulon$target[regulon$tf=="A"][c(1,2,4)]
regulon$target[regulon$tf=="C"][6] <- regulon$target[regulon$tf=="A"][3]
regulon$target[is.na(regulon$target)] <- sample(LETTERS[10:14], sum(is.na(regulon$target)), replace = TRUE)
regulon$weights <- runif(nrow(regulon))

res <- vector(mode = "list", 3)
names(res) <- LETTERS[2:4]
res <- lapply(res, function(x) data.frame(target=character(0), focal_weight=numeric(0),
                                          other_tf_weight=numeric(0), weight_product=numeric(0)))

res$B <- rbind(res$B, data.frame(target=regulon$target[regulon$tf=="A"][c(1,2,4)],
                                 focal_weight=regulon$weights[regulon$tf=="A"][c(1,2,4)],
                                 other_tf_weight=regulon$weights[regulon$tf=="B"][c(1,3,5)],
                                 weight_product=NA))

res$C <- rbind(res$C, data.frame(target=regulon$target[regulon$tf=="A"][3],
                                 focal_weight=regulon$weights[regulon$tf=="A"][3],
                                 other_tf_weight=regulon$weights[regulon$tf=="C"][6],
                                 weight_product=NA))

res <- lapply(res, function(x) {x$weight_product <- x$focal_weight*x$other_tf_weight; x})

test_graph <- buildGraph(regulon)
partners_graph <- findPartners(test_graph, "A")

test_that("findPartners function works correctly", {
  expect_identical(partners_graph, res)
})
