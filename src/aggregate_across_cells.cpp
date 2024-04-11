#include "Rtatami.h"
#include "AggregateAcrossCells.h"

//[[Rcpp::export(rng=false)]]
SEXP aggregate_across_cells(SEXP x, Rcpp::List groupings, int nthreads) {
  Rtatami::BoundNumericPointer converted(x);
  const auto& mat = converted->ptr;
  size_t NC = mat->ncol();
  
  size_t nfactors = groupings.size();
  std::vector<Rcpp::IntegerVector> groups(nfactors);
  if (nfactors == 0) {
    throw std::runtime_error("expected a non-zero number of grouping factors");
  }
  
  std::vector<const int*> gptrs(nfactors);
  for (size_t g = 0; g < nfactors; ++g) {
    groups[g] = Rcpp::IntegerVector(groupings[g]);
    if (groups[g].size() != NC) {
      throw std::runtime_error("each factor should have the same length as the number of columns");
    }
    gptrs[g] = static_cast<const int*>(groups[g].begin());
  }
  
  Rcpp::IntegerVector fac(NC);
  int * fptr = fac.begin();
  auto combined = scran::AggregateAcrossCells::combine_factors(NC, gptrs, fptr);
  
  // Copying things over.
  Rcpp::IntegerVector counts(combined.counts.begin(), combined.counts.end());
  Rcpp::List factors(nfactors);
  for (size_t f = 0; f < nfactors; ++f) {
    const auto& current = combined.factors[f];
    factors[f] = Rcpp::IntegerVector(current.begin(), current.end());
  }
  
  // Setting up the output matrices.
  size_t NR = mat->nrow();
  size_t ncombos = combined.counts.size();
  Rcpp::NumericMatrix sums(NR, ncombos);
  Rcpp::IntegerMatrix detected(NR, ncombos);
  
  std::vector<double*> sptrs;
  std::vector<int*> dptrs;
  for (size_t i = 0; i < ncombos; ++i) {
    sptrs.push_back(static_cast<double*>(sums.begin() + i * NR));
    dptrs.push_back(static_cast<int*>(detected.begin() + i * NR));
  }
  
  scran::AggregateAcrossCells runner;
  runner.set_num_threads(nthreads);
  runner.run(mat.get(), fptr, sptrs, dptrs);
  
  return Rcpp::List::create(
    Rcpp::Named("sums") = sums,
    Rcpp::Named("detected") = detected,
    Rcpp::Named("combinations") = factors,
    Rcpp::Named("counts") = counts,
    Rcpp::Named("index") = fac
  );
}