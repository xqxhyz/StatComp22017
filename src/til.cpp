#include <Rcpp.h>
using namespace Rcpp;

//' @title Compute TIL Score using Rcpp
//' @param LYM the number of patches of LYM
//' @param STR the number of patches of STR
//' @param TUM the number of patches of TUM
//' @param TIL3 the TIL Score of LYM
//' @param TIL7 the TIL Score of STR
//' @param TIL8 the TIL Score of TUM
//' @useDynLib StatComp22017
//' @export
// [[Rcpp::export]]
double til(double LYM, double STR, double TUM, double TIL3, double TIL7, double TIL8){
  double k3 = LYM / (STR+TUM+LYM);
  double k7 = STR / (STR+TUM+LYM);
  double k8 = TUM / (STR+TUM+LYM);
  double til_t = k3*TIL3+k7*TIL7+k8*TIL8;
  return til_t;
}
