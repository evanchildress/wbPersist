#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector grow(NumericMatrix pars,
                       IntegerVector river,
                       NumericVector startLength,
                       NumericVector flow,
                       NumericVector biomass,
                       NumericVector performance) {

  int n = river.size();
  NumericVector endLength(n);
  for(int i = 0; i<n; ++i) {
    endLength[i] = (pars(1,river[i]) +
                      pars(2,river[i])*startLength[i] +
                      pars(3,river[i])*flow[river[i]] +
                      pars(4,river[i])*biomass[river[i]])*
                   performance[i];
  }
  return endLength;
}
