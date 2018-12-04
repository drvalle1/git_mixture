#include <Rcpp.h>
#include <iostream>
#include <ctime>
#include <fstream>
using namespace Rcpp;

/***************************************************************************************************************************/
/*********************************                      UTILS          *****************************************************/
/***************************************************************************************************************************/

//' Identify for which index the cumulative sum of prob is smaller than value
//' 
//' Identify for which index the cumulative sum of probabilities (prob) is smaller than
//' a given value (value). This function helps with multinomial draws
//' 
//' @param value a real number between 0 and 1
//' @param prob a vector of probabilities that sum to one
//' @return this function returns the integer res
//' @export

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

//' Generate samples from a multinomial distribution with n=1
//' 
//' Generates samples from a multinomial distribution with n=1. 
//' Number of samples are equal to the number of rows in the matrix prob
//' Number of classes are equal to the number of columns in the matrix prob
//' 
//' @param prob L x K matrix containing the probability for each location l and class k
//'        Probabilities sum to one across rows
//' @param randu set of uniform random variables
//' @return this function returns a vector of length L, where cs[i] contains
//'         a random variable from a multinomial distribution with n=1 and probability = prob(i,_)
//' @export
// [[Rcpp::export]]
IntegerVector rmultinom1(NumericMatrix prob, NumericVector randu) {
  
  IntegerVector cs(prob.nrow());

  for(int i=0; i<prob.nrow();i++){
    cs[i]=whichLessDVPresence(randu[i],prob(i,_));
  }
  return cs;
}

//' Summarize the data
//' 
//' This function summarizes the data by calculating the number of locations in each 
//' group and species for two cases: dat(i,j)=1 (stored in res1) and dat(i,j)=0 (stored in res0)
//' 
//' @param dat this matrix has L rows (e.g., locations) and S columns (e.g., species)
//'        and contains the presence-absence data
//' @param z vector with cluster assignment for each location
//' @param nspp number of species
//' @param nloc number of locations
//' @param ngroup maximum number of groups
//' @return this function returns a list containing two matrices: res1 and res0
//' @export

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
