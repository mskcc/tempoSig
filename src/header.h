struct Param{
  const Rcpp::NumericVector &D;
  const Rcpp::NumericMatrix &ref;
};

double lnl(const gsl_vector *v, void *params);

void dlnl(const gsl_vector *v, void *params, gsl_vector *df);

void lnl_dlnl(const gsl_vector *v, void *params, double *f, gsl_vector *df);
