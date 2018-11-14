#include <Rcpp.h>
#include <iostream>
#include <ctime>
#include <fstream>
using namespace Rcpp;

/***************************************************************************************************************************/
/*********************************                      UTILS          *****************************************************/
/***************************************************************************************************************************/

// This function helps with multinomial draws
int whichLessDVPresence(double value, NumericVector prob) {
  int res=-1;
  double probcum = 0;
  
  for (int i = 0; i < prob.length(); i++) {
    probcum = probcum + prob(i);
    if (value < probcum) {
      res = i;
      break;
    }
  }
  return res;
}

//' This function samples cs's
// [[Rcpp::export]]
IntegerVector rmultinom1(NumericMatrix prob, NumericVector randu) {
  
  IntegerVector cs(prob.nrow());

  for(int i=0; i<prob.nrow();i++){
    cs[i]=whichLessDVPresence(randu[i],prob(i,_));
  }
  return cs;
}

//' This function calculates ncs1 and ncs0
// [[Rcpp::export]]
Rcpp::List ncs(IntegerMatrix dat, IntegerVector z, 
               int nspp, int nloc, int ngroup) {
  
  IntegerMatrix res1(ngroup,nspp);
  IntegerMatrix res0(ngroup,nspp);
  
  for(int i=0; i<nloc; i++){
    for(int j=0; j<nspp; j++){
      if(dat(i,j)==1) res1(z(i),j)=res1(z(i),j)+1;
      if(dat(i,j)==0) res0(z(i),j)=res0(z(i),j)+1;
    }
  }

  Rcpp::List resTemp = Rcpp::List::create(Rcpp::Named("ncs1") = res1,
                                          Rcpp::Named("ncs0") = res0);
  return(resTemp);
}
