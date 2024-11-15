addPadding <- function(gr){ # extend regions covered by chip-seq to generate partially overlapping peaks
    if(length(gr)<2) return(resize(gr, width=width(gr)+1000, fix = "center"))
    gr <- reduce(gr)
    distances_downstream <- c(500, distance(gr[1:((length(gr)-1))], gr[2:length(gr)]))
    distances_upstream <- c(distances_downstream[2:length(distances_downstream)], 500)
    distances_downstream <- sapply(floor(distances_downstream/2), min, 500) # avoid merging adjacent binding regions during subsequent operations
    distances_upstream <- sapply(floor(distances_upstream/2), min, 500)
    distances_downstream[is.na(distances_downstream)] <- 500
    distances_upstream[is.na(distances_upstream)] <- 500
    gr <- resize(gr, width=width(gr)+distances_downstream, fix = "end")
    gr <- resize(gr, width=width(gr)+distances_upstream, fix = "start")
    gr
}

create_peak_sequence <- function(gr){ # narrow range to the maximum size of 501 bp
    start_positions <- width(gr)-501
    if(start_positions<1) return(gr)
    else {
        start_pos <- start(gr)+sample(0:start_positions,1)
        end_pos <- start_pos+500
    }
    return(GRanges(ranges=IRanges::IRanges(start = start_pos, end=end_pos, seqnames=seqnames(gr))))
}

select_hits <- function(o){ # choose random region if multiple regions overlap with the same TF binding site
    if(length(o)==0) return(numeric(0))
    first_occurrence <- which(!duplicated(subjectHits(o)))
    last_occurrence <- (length(o)+1) - which(!duplicated(rev(subjectHits(o))))
    ind <- first_occurrence + sapply(last_occurrence - first_occurrence, function(x) sample(0:x,1))
    queryHits(o)[ind]
}

set.seed(3021)
grl_binding <- getTFMotifInfo()
grl_binding <- grl_binding[sample(1:length(grl_binding), 30)]
grl_binding <- GRangesList(lapply(grl_binding, function(x) x[sample(length(x), round(0.2*length(x)))]))
grl_binding_copy <- grl_binding
selected_tfs <- sort(sample(length(grl_binding), round(0.8*length(grl_binding)))) # choose TFs to be covered by peaks
grl_binding_negative <- GRangesList(lapply(grl_binding, function(x) gaps(x)))
unoccupied_sequences <- unlist(GRangesList(S4Vectors::Reduce(GenomicRanges::intersect, grl_binding_negative)))
pairs <- matrix(sample(selected_tfs, 18), ncol=2)

shared_regions <- list()

# look for the regions covered by selected pair of TFs and not covered by any other TF bindinf site

for(i in 1:nrow(pairs)){
    unoccupied_by_other_tfs <- S4Vectors::Reduce(GenomicRanges::intersect, grl_binding_negative[-pairs[i,]])
    shared_regions[[i]] <- GenomicRanges::intersect(grl_binding[[pairs[i,1]]], GenomicRanges::intersect(grl_binding[[pairs[i,2]]], unoccupied_by_other_tfs))
}

shared_regions <- GRangesList(shared_regions)

overlaps <- data.frame()
peaks <- GRanges()

for(i in seq_along(shared_regions)){
    n_peaks <- min(length(shared_regions[[i]]), 4)
    if(n_peaks>0){
        regions <- sample(length(shared_regions[[i]]), n_peaks)
        for(region in regions){
            peaks<-c(peaks, create_peak_sequence(shared_regions[[i]][region]))
            overlaps <- rbind(overlaps, data.frame(idxATAC = length(peaks), idxTF = pairs[i,], tf = names(grl_binding)[pairs[i,]]))
        }
    }
}

# look for regions covered only by selected TFs

grl_binding_unique <- GRangesList()
unoccupied_by_other_tfs <- vector(mode = "list", length = length(grl_binding))
for(i in selected_tfs){
    if(is.null(unoccupied_by_other_tfs[[i]]))  unoccupied_by_other_tfs[[i]] <- S4Vectors::Reduce(GenomicRanges::intersect, grl_binding_negative[-i])
}
# clean by subtracting sequences covered by other TFs
grl_binding_unique <- GRangesList(lapply(selected_tfs, function(i,grl, grl_uo) GenomicRanges::intersect(grl[[i]], grl_uo[[i]]), grl=grl_binding, grl_uo=unoccupied_by_other_tfs))
grl_binding_extended <- GRangesList(lapply(grl_binding_unique, addPadding))
grl_binding_extended <- GRangesList(lapply(grl_binding_extended, function(x) GenomicRanges::intersect(x, unoccupied_sequences))) # extension should not overlap with bindings sites of other TFs
# glue paddings with cleaned TF binding sites
grl_binding_extended <- GRangesList(lapply(selected_tfs, function(i, grl,grl_ex) reduce(c(grl[[which(selected_tfs==i)]], grl_ex[[which(selected_tfs==i)]])),
                                      grl = grl_binding_unique, grl_ex = grl_binding_extended))
grl_binding_extended <- GRangesList(lapply(selected_tfs, function(i, grl,grl_ex) {o <- findOverlaps(grl_ex[[which(selected_tfs==i)]], grl[[i]]);grl_ex[[which(selected_tfs==i)]][select_hits(o)]},
                                      grl = grl_binding, grl_ex = grl_binding_extended))

# save extended binding sites as peaks and register them in overlaps object
new_peaks <- GRanges()
for(tf_id in seq_along(grl_binding_extended)){
    regions_n <- length(grl_binding_extended[[tf_id]])
    new_peaks_n <- min(5, sample(regions_n, 1)) # select number of peaks to overlap with the current TF binding sites
    if(new_peaks_n > 0){
        original_TF_ID <- selected_tfs[tf_id]
        overlaps <- rbind(overlaps, data.frame(idxATAC = length(peaks)+length(new_peaks)+seq_len(new_peaks_n),
                                               idxTF = original_TF_ID, tf = names(grl_binding)[original_TF_ID]))
        new_peaks <- c(new_peaks, grl_binding_extended[[tf_id]][sample(regions_n,new_peaks_n)])
    }
}

to_narrow <- which(width(new_peaks)>501)
if(length(to_narrow)>0){
    new_start <- sapply(width(new_peaks[to_narrow])-501, function(x) sample(0:x,1)) + start(new_peaks[to_narrow])
    start(new_peaks)[to_narrow] <- new_start
    end(new_peaks)[to_narrow] <- new_start+500
}

peaks <- c(peaks, new_peaks)

# add peaks with no overlaps
peaks <- c(peaks, unoccupied_sequences[sample(length(unoccupied_sequences), 30)])
# shuffle peak order
new_order <- sample(length(peaks))
peaks <- peaks[new_order]
overlaps$idxATAC <- match(overlaps$idxATAC, new_order)
p2g <- data.frame(idxATAC = overlaps$idxATAC)
peakMatrix <- SingleCellExperiment(assays = list(counts = matrix(sample(c(0,1), length(peaks)*200,replace = TRUE),length(peaks),200)), rowRanges = peaks)
overlaps <- overlaps[order(overlaps$idxATAC, overlaps$idxTF),]
rownames(overlaps) <- 1:nrow(overlaps)

overlaps_test <- addTFMotifInfo(p2g, grl=grl_binding, peakMatrix = peakMatrix)

test_that("addTFMotifInfo works correclty", {
    expect_identical(overlaps_test, overlaps)
})
#################

grl <- getTFMotifInfo()

test_that("getTFMotifInfo works with default settings", {
    expect_equal(length(grl), 1558)
    expect_equal(length(unlist(grl)), 54501169)
    expect_equal(sum(width(unlist(grl))), 21165342299)
})


grl <- getTFMotifInfo(genome = "hg38", source = "atlas.sample")

test_that("getTFMotifInfo works for chip-atlas sample specific data", {
    expect_equal(length(grl), 1019)
    expect_equal(length(unlist(grl[[1]])), 43050)
    expect_equal(sum(width(unlist(grl[[8]]))), 6644003)
})

set.seed(4722)
all_peaks <- 1:1000
all_genes <- paste0("gene_", 1:3000)
peak_gene_links <- data.frame(idxATAC = sample(all_peaks, 4000, replace = TRUE), target = sample(all_genes, 4000, replace = TRUE))
peak_gene_links <- peak_gene_links[!duplicated(peak_gene_links),]
peak_gene_links$Correlation.all <- 2*runif(nrow(peak_gene_links))-1
peak_gene_links <- S4Vectors::DataFrame(peak_gene_links)
tf_to_region_links <- data.frame(TF = c(sample(paste0("gene_", 1:10), 60, replace = TRUE),
                                        sample(paste0("tf_", 1:50), 1000, replace = TRUE)),
                                 idxATAC = sample(1:2000, 1060, replace = TRUE))
res <- S4Vectors::merge(peak_gene_links, tf_to_region_links, by = "idxATAC")
corr_matrix <- as.matrix(res[,"Correlation.all", drop = FALSE])
colnames(corr_matrix) <- "all"
res$corr <-  corr_matrix
res[["Correlation.all"]] <- NULL
regulon <- getRegulon(peak_gene_links, tf_to_region_links)

test_that("getRegulon works correctly", {
    expect_identical(regulon, res)
})
