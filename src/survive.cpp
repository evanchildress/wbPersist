#include <Rcpp.h>
using namespace Rcpp;

//' Gets a survival probability given the environmental contribution, start length, and parameters
//'
//' @param river Vector of river locations of individuals
//' @param stage Vector of stages
//' @param forkLength Vector of lengths of individuals
//' @param phiBeta5 Matrix of length betas
//' @param envLogitPhi Matrix of environmental contribution to surival from \code{getEnvLogitPhi()}
//' @export
//' @useDynLib wbPersist
//' @importFrom Rcpp sourceCpp
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
