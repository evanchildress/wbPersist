#include <Rcpp.h>
using namespace Rcpp;

//'Estimates actual growth given performance, other covariates, and parameters
//'
//'@param pars Matrix of betas from growth model
//'@param river Vector of river locations for individuals
//'@param startLength Vector of starting lengths
//'@param flow Vector of mean flow for the growth period
//'@param biomass Vector of brook trout biomass during the growth period
//'@param performance Vector of growth performance from the \code{getGrowthPerformance()}
//'@export
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
    endLength[i] = startLength[i] +
                   (pars(0,river[i]-1) +
                      pars(1,river[i]-1)*startLength[i] +
                      pars(2,river[i]-1)*flow[river[i]-1] +
                      pars(3,river[i]-1)*biomass[river[i]-1])*
                   performance[river[i]-1];
  }
  return endLength;
}
