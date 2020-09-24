#include <Rcpp.h>
#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>
#include "header.h"

using namespace Rcpp;

// [[Rcpp::export]]
List mlestimate(NumericVector D, NumericVector x0, NumericMatrix ref, IntegerVector Itmax,
                NumericVector Tol){

  int N = x0.size();  // no. of signatures

  unsigned int Imax = Itmax[0];

  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s;
  gsl_vector *x;
  gsl_multimin_function_fdf func;

  func.n = N;
  func.f = lnl;
  func.df = dlnl;
  func.fdf = lnl_dlnl;

  x = gsl_vector_alloc(N);
  T = gsl_multimin_fdfminimizer_vector_bfgs2;
  s = gsl_multimin_fdfminimizer_alloc(T, N);

  Param par = {D, ref};
  func.params = &par;
  for(int i=0; i<N; i++)
    gsl_vector_set(x, i, x0[i]); // initial guess

  gsl_multimin_fdfminimizer_set(s, &func, x, 0.1, 0.1);

  size_t iter=0;
  int status;
  do{
    iter++;
    status =  gsl_multimin_fdfminimizer_iterate(s);
    if(status){
      Rcpp::Rcerr << " GSL status code " << status << std::endl;
      break;
    }
    status = gsl_multimin_test_gradient(s->gradient, Tol[0]);
  }while(status == GSL_CONTINUE && iter < Imax);
  if(iter == Imax)
    Rcpp::Rcerr << "BFGS2 iteration maximum reached\n";

  std::vector<double> v(N);
  for(int i=0; i<N; i++)
    v[i] = gsl_vector_get(s->x, i);

  List f = List::create(Named("x") = v, Named("lkh") = -(s->f));

  gsl_multimin_fdfminimizer_free(s);
  gsl_vector_free(x);

  return f;
}

double lnl(const gsl_vector *v, void *params){

  Param *par = (Param *)params;
  NumericVector D = par->D;
  NumericMatrix sig = par->ref;

  int M = sig.nrow();
  int N = sig.ncol();

  std::vector<double> x(N);
  double p0 = 0.0;
  for(int k=0; k<N; k++){
    x[k] = gsl_vector_get(v, k);
    p0 += x[k]*x[k];
  }

  double E = 0.0;
  double sum = 0.0;

  for(int i=0; i<M-1; i++){
    double p1 = 0;
    for(int k=0; k<N; k++)
      p1 += sig(i, k)*x[k]*x[k];
    E -= D[i]*log(p1/p0);
    sum += p1/p0;
  }
//if(sum > 1.0 || sum < 0){
  if(sum > 1.000001 || sum < 0)
    Rcpp::stop("Error in lnl\n");
  if(sum < 1.0)
    E -= D[M-1]*log(1-sum);
  else E -= D[M-1]*log(sum-1);

  return E;
}

void dlnl(const gsl_vector *v, void *params, gsl_vector *df){

  Param *par = (Param *)params;
  NumericVector D = par->D;
  NumericMatrix sig = par->ref;

  int M = sig.nrow();
  int N = sig.ncol();

  std::vector<double> x(N);
  double p0 = 0.0;
  for(int k=0; k<N; k++){
    x[k] = gsl_vector_get(v, k);
    p0 += x[k]*x[k];
  }
  std::vector<double> p1(M-1);
  double pL = 0.0;
  for(int i=0; i<M-1; i++){
    p1[i] = 0.0;
    for(int k=0; k<N; k++)
      p1[i] += sig(i,k)*x[k]*x[k];
    p1[i] /= p0;
    pL += p1[i];
  }
  pL = 1.0 - pL;

  for(int k=0; k<N; k++){
    double sm = 0.0;
    for(int i=0; i<M-1; i++)
      sm += (-D[i]/p1[i] + D[M-1]/pL)*(sig(i,k)-p1[i]);
    sm *= 2.0*x[k]/p0;
    gsl_vector_set(df, k, sm);
  }
}

void lnl_dlnl(const gsl_vector *v, void *params, double *f, gsl_vector *df){

  *f = lnl(v, params);
  dlnl(v, params, df);

}
