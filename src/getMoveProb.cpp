#include <Rcpp.h>
using namespace Rcpp;


//' Return the probability of moving to each river given a starting location
//'
//' Requires that movement probabilities are given as a nRivers x nRivers matrix, so this needs to be run for a single season.
//'
//' @param river Vector of starting locations
//' @param moveProb Matrix of season-specific movement probabilities
//' @export
// [[Rcpp::export]]
NumericMatrix getMoveProb(NumericVector river,
                   NumericMatrix moveProb) {

  int n = river.size();
  NumericMatrix m(n,4);

  for(int i = 0; i<n; ++i){
    m(i,_) = moveProb(river[i]-1,_);
  }

  return m;
}
