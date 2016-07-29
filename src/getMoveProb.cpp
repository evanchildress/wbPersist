#include <Rcpp.h>

using namespace Rcpp;

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
