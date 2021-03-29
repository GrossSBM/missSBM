#ifndef _utils_optim_H
#define _utils_optim_H

inline double logistic(double x) {
  return (1/(1+std::exp(-x)));
}

inline arma::mat logistic(arma::mat & x) {
  return (1/(1+exp(-x)));
}

inline double g(double x) {
  return (-std::log(1+std::exp(-x)));
}

inline double g_prime(double x) {
  return (1/(1+std::exp(x)));
}

inline double u(double x) {
  return (1/(1+std::exp(x)));
}

inline arma::vec softmax(arma::vec &x) {
  double b = max(x) ;
  return(exp(x - b) / sum(exp(x - b))) ;
}

#endif
