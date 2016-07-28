#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector survive(NumericVector river,
                      NumericVector stage,
                      NumericVector forkLength,
                      NumericMatrix phiBeta5,
                      NumericMatrix envLogitPhi) {
  
  int n = river.size();
  NumericVector logitPhi(n);
  NumericVector phi(n);
  
  for(int i = 0; i<n; ++i){
    logitPhi[i] = envLogitPhi(river[i],stage[i]) + 
                  phiBeta5(river[i],stage[i])*forkLength[i];
  }
  
  phi = 1/(1+exp(-logitPhi));
  return phi;
}
