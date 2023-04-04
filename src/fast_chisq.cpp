#include "Rcpp.h"

//[[Rcpp::export(rng=false)]]
Rcpp::List fast_chisq(
    Rcpp::IntegerVector peak_ordered, 
    Rcpp::IntegerVector tf_by_peak,
    Rcpp::IntegerVector target_by_peak, 
    Rcpp::IntegerVector target_ordered,

    int npeaks,
    Rcpp::NumericVector peakmat_x, // tatami in flux so we'll just drag the vectors in.
    Rcpp::IntegerVector peakmat_i, 
    Rcpp::IntegerVector peakmat_p, 
    double peak_cutoff,

    int ngenes,
    Rcpp::NumericVector expmat_x,
    Rcpp::IntegerVector expmat_i, 
    Rcpp::IntegerVector expmat_p,
    double exp_cutoff,

    int nclusters,
    Rcpp::IntegerVector clusters)
{
    size_t nrows = peak_ordered.size();
    if (nrows != tf_by_peak.size()) {
        throw std::runtime_error("'peak_ordered' and 'tf_by_peak' should have the same length");
    }
    if (nrows != target_by_peak.size()) {
        throw std::runtime_error("'peak_ordered' and 'target_by_peak' should have the same length");
    }
    if (nrows != target_ordered.size()) {
        throw std::runtime_error("'peak_ordered' and 'target_ordered' should have the same length");
    }

    // Building a reverse index on the peaks and targets.
    std::vector<size_t> peak_start(npeaks), peak_end(npeaks);
    for (size_t i = 0; i < nrows; ++i) {
        auto& current_start = peak_start[peak_ordered[i]];
        auto& current_end = peak_end[peak_ordered[i]];
        if (current_end == 0) {
            current_start = i;
            current_end = i + 1;
        } else {
            ++current_end;
        }

        if (i && peak_ordered[i] < peak_ordered[i-1]) {
            throw std::runtime_error("'peak_ordered' should be in ascending order");
        }
    }

    std::vector<size_t> target_start(ngenes), target_end(ngenes);
    for (size_t i = 0; i < nrows; ++i) {
        auto& current_start = target_start[target_ordered[i]];
        auto& current_end = target_end[target_ordered[i]];
        if (current_end == 0) {
            current_start = i;
            current_end = i + 1;
        } else {
            ++current_end;
        }

        if (i && target_ordered[i] < target_ordered[i-1]) {
            throw std::runtime_error("'target_ordered' should be in ascending order");
        }
    }

    // Checking the number of columns.
    int ncols = expmat_p.size() - 1;
    if (expmat_p.size() != peakmat_p.size()) {
        throw std::runtime_error("expession and peak matrices should have the same number of columns");
    }

    std::vector<unsigned char> exists_in_exp(ngenes);
    Rcpp::NumericMatrix output_triple(nrows, nclusters);
    Rcpp::NumericMatrix output_peak(nrows, nclusters);
    Rcpp::NumericMatrix output_target(nrows, nclusters);

    for (int c = 0; c < ncols; ++c) {
        int clust = clusters[c];
        double* trptr = output_triple.begin() + clust * nrows;
        double* pkptr = output_peak.begin() + clust * nrows;
        double* taptr = output_target.begin() + clust * nrows;

        // First pass to tag each gene for whether it's expressed.
        int exstart = expmat_p[c];
        int exend = expmat_p[c + 1];
        for (int i = exstart; i < exend; ++i) {
            if (expmat_x[i] > exp_cutoff) {
                int index = expmat_i[i];
                exists_in_exp[index] = 1;
            }
        }

        // Pass through the peak matrix and collect the peak + TF and peak + TF + target counts.
        int pkstart = peakmat_p[c];
        int pkend = peakmat_p[c + 1];
        for (int i = pkstart; i < pkend; ++i) {
            if (peakmat_x[i] > peak_cutoff) {
                int index = peakmat_i[i];
                auto regstart = peak_start[index];
                auto regend = peak_end[index];

                for (size_t j = regstart; j < regend; ++j) {
                    if (exists_in_exp[tf_by_peak[j]]) {
                        ++pkptr[j];
                        if (exists_in_exp[target_by_peak[j]]) {
                            ++trptr[j];
                        }
                    }
                }
            }
        }

        // Pass through the expression matrix and collect target counts.
        for (int i = exstart; i < exend; ++i) {
            if (expmat_x[i] > exp_cutoff) {
                int index = expmat_i[i];
                auto regstart = target_start[index];
                auto regend = target_end[index];

                for (size_t j = regstart; j < regend; ++j) {
                    ++taptr[j];
                }

                // Also resetting the buffer while we're here.
                exists_in_exp[index] = 0;
            }
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("triple") = output_triple,
        Rcpp::Named("peak") = output_peak,
        Rcpp::Named("target") = output_target
    );
}
