#include <Rcpp.h>
#include <iostream>
#include <ctime>
#include <fstream>
using namespace Rcpp;

/***************************************************************************************************************************/
/*********************************                      UTILS          *****************************************************/
/***************************************************************************************************************************/

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
