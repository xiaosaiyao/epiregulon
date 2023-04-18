#include "Rcpp.h"

#include <algorithm>
#include <vector>

struct ComputeWorkspace {
    ComputeWorkspace(size_t n) : less_than(n), equal0(n), equal1(n) {}
    std::vector<double> less_than;
    std::vector<double> equal0;
    std::vector<double> equal1;
};

void compute_u_statistic(
    const std::vector<std::pair<double, int> >& input, 
    const std::vector<double>& num_zeros0, 
    const std::vector<double>& num_zeros1,
    const int* clusters,
    const unsigned char* okay,
    ComputeWorkspace& work,
    double* output)
{
    int nclusters = num_zeros0.size();
    auto& less_than = work.less_than;
    std::fill(less_than.begin(), less_than.end(), 0);

    // Values == 0.
    for (int c = 0; c < nclusters; ++c) {
        if (num_zeros1[c]) {
            output[c] += num_zeros1[c] * (less_than[c] + 0.5 * num_zeros0[c]);
        }
        less_than[c] += num_zeros0[c];
    }

    // Values > 0.
    size_t pos = 0;
    auto& equal0 = work.equal0;
    std::fill(equal0.begin(), equal0.end(), 0);
    auto& equal1 = work.equal1;
    std::fill(equal1.begin(), equal1.end(), 0);

    while (pos != input.size()) {
        const auto& current = input[pos];

        ++pos;
        bool tied = false;
        while (pos != input.size() && input[pos].first == current.first) {
            tied = true;
            int c = input[pos].second;
            auto& equal = (okay[c] ? equal1 : equal0);
            ++equal[clusters[c]];
            ++pos;
        }

        int self = current.second;
        int self_okay = okay[self];
        int self_clust = clusters[self];

        if (tied) {
            // Adding the current element to the tie.
            auto& equal = (self_okay ? equal1 : equal0);
            ++equal[self_clust]; 

            for (int c = 0; c < nclusters; ++c) {
                if (equal1[c]) {
                    output[c] += equal1[c] * (less_than[c] + 0.5 * equal0[c]);
                }
                less_than[c] += equal0[c];
            }

            std::fill(equal0.begin(), equal0.end(), 0);
            std::fill(equal1.begin(), equal1.end(), 0);

        } else if (self_okay) {
            output[self_clust] += less_than[self_clust];
            ++less_than[self_clust];
        }
    };

    return;
}

//[[Rcpp::export(rng=false)]]
Rcpp::List fast_wilcox(
    Rcpp::NumericVector exprs_x,
    Rcpp::IntegerVector exprs_i,
    Rcpp::IntegerVector exprs_p,
    Rcpp::NumericVector peak_x,
    Rcpp::IntegerVector peak_i,
    Rcpp::IntegerVector peak_p,
    Rcpp::IntegerVector target_id,
    Rcpp::IntegerVector tf_id,
    Rcpp::IntegerVector peak_id,
    Rcpp::IntegerVector clusters,
    int ngroups)
{
    size_t ncells = clusters.size();

    int last = -1;
    std::vector<std::pair<double, int> > sortspace;
    sortspace.reserve(ncells);
    std::vector<std::vector<std::pair<double, int> > > cluster_space(ngroups);

    std::vector<unsigned char> okay(ncells);
    std::vector<int> okay_indices;

    std::vector<double> okay_zeros(ngroups), okay_total(ngroups), notokay_zeros(ngroups), notokay_total(ngroups);
    ComputeWorkspace workspace(ngroups);

    std::vector<int> full_cluster(ncells);
    std::vector<double> full_okay_zeros(1), full_notokay_zeros(1);
    ComputeWorkspace full_workspace(1);

    size_t nregs = target_id.size();
    Rcpp::NumericMatrix output_u(ngroups + 1, nregs);
    Rcpp::NumericMatrix output_t0(ngroups + 1, nregs);
    Rcpp::NumericMatrix output_t1(ngroups + 1, nregs);

    for (size_t i = 0; i < nregs; ++i) {
        // Only resorting if we've moved onto a different target gene.
        if (last != target_id[i]) {
            sortspace.clear();
            int offset = exprs_p[target_id[i]], last = exprs_p[target_id[i]+1];
            for (int j = offset; j < last; ++j) {
                if (exprs_x[j] > 0) {
                    sortspace.emplace_back(exprs_x[j], exprs_i[j]);
                }
            }
            std::sort(sortspace.begin(), sortspace.end());
            last = target_id[i];
        }

        // Identifying the cells with TF expression + peak.
        {
            int tf_offset = exprs_p[tf_id[i]], tf_last = exprs_p[tf_id[i]+1];
            int peak_offset = exprs_p[peak_id[i]], peak_last = exprs_p[peak_id[i]+1];
            int k = peak_offset;

            for (int j = tf_offset; j < tf_last; ++j) {
                if (exprs_x[j] <= 0) {
                    continue;
                }

                int icurrent = exprs_i[j];
                while (k < peak_last && peak_i[k] < icurrent) {
                    ++k;
                }

                if (k == peak_last) {
                    break;
                }

                if (peak_i[k] == icurrent && peak_x[k] > 0) {
                    okay_indices.push_back(icurrent);
                    okay[icurrent] = 1;
                }
            }
        }

        {
            std::fill(okay_zeros.begin(), okay_zeros.end(), 0);
            std::fill(okay_total.begin(), okay_total.end(), 0);
            std::fill(notokay_zeros.begin(), notokay_zeros.end(), 0);
            std::fill(notokay_total.begin(), notokay_total.end(), 0);

            // Looping through the sortspace and creating cluster-specific sorted vectors.
            for (const auto& ss : sortspace) {
                int c = clusters[ss.second];
                if (okay[ss.second]) {
                    ++okay_zeros[c];
                    ++okay_total[c];
                } else {
                    ++notokay_zeros[c];
                    ++notokay_total[c];
                }
            }

            // Computing the U statistic between the cells with TF+peak, and those without both.
            for (int c = 0; c < ngroups; ++c) {
                okay_zeros[c] = okay_total[c] - okay_zeros[c];
                notokay_zeros[c] = notokay_total[c] - notokay_zeros[c];
            }

            auto col_u = output_u.column(i);
            auto output_u_ptr = static_cast<double*>(col_u.begin());
            compute_u_statistic(sortspace, notokay_zeros, okay_zeros, clusters.begin(), okay.data(), workspace, output_u_ptr);

            full_okay_zeros[0] = std::accumulate(okay_zeros.begin(), okay_zeros.end(), 0);
            full_notokay_zeros[0] = std::accumulate(notokay_zeros.begin(), notokay_zeros.end(), 0);
            compute_u_statistic(sortspace, full_okay_zeros, full_notokay_zeros, full_cluster.data(), okay.data(), full_workspace, output_u_ptr + ngroups);

            // Copying the totals to the output.
            auto col_t0 = output_t0.column(i);
            std::copy(notokay_total.begin(), notokay_total.end(), col_t0.begin());
            col_t0[ngroups] = std::accumulate(okay_total.begin(), okay_total.end(), 0);

            auto col_t1 = output_t1.column(i);
            std::copy(okay_total.begin(), okay_total.end(), col_t1.begin());
            col_t1[ngroups] = std::accumulate(notokay_total.begin(), notokay_total.end(), 0);
        }

        // Mopping up for the next regulon triplet.
        {
            for (auto o : okay_indices) {
                okay[o] = 0;
            }
            okay_indices.clear();

            for (auto& cs : cluster_space) {
                cs.clear();
            }
        }

    }

    return Rcpp::List::create(
        Rcpp::Named("U") = output_u,
        Rcpp::Named("T0") = output_t0,
        Rcpp::Named("T1") = output_t1
    );
}
