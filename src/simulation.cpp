#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


// [[Rcpp::export(.garchsimvec)]]
List garchsimvec(arma::mat& epsilon, arma::mat& sigma_sqr_sim, const arma::mat& z,
                 const arma::vec& variance_intercept, const arma::mat& init, const arma::vec& alpha,
                 const arma::vec& beta, const double mu, const arma::ivec& order) {
    const int maxpq = (int) max(order);
    int h = z.n_cols - maxpq; // Assuming z already includes space for burn-in
    int nsim = z.n_rows;
    arma::mat sigma_sim(nsim, h + maxpq);
    arma::mat series_sim(nsim, h + maxpq);
    int i,j = 0;
    for(i=maxpq;i<(maxpq+h);i++) {
        sigma_sqr_sim.col(i).fill(arma::as_scalar(variance_intercept(i)));
        if (order(0) > 0) {
            for(j=0;j<order(0);j++) {
                if((order(0) + j) >= i) {
                    sigma_sqr_sim.col(i) += alpha(j) * init.col(j);
                } else {
                    sigma_sqr_sim.col(i) += alpha(j) * pow(epsilon.col(i - j - 1), 2.0);
                }
            }
        }
        if (order(1) > 0) {
            for(j=0;j<order(1);j++) {
                sigma_sqr_sim.col(i) += beta(j) * sigma_sqr_sim.col(i - j - 1);
            }
        }
        sigma_sim.col(i) = sqrt(sigma_sqr_sim.col(i));
        epsilon.col(i) = z.col(i) % sigma_sim.col(i);
        series_sim.col(i) = epsilon.col(i) + mu;
    }
    return List::create(Named("sigma") = sigma_sim, Named("series") = series_sim, Named("epsilon") = epsilon);
}

// [[Rcpp::export(.egarchsimvec)]]
List egarchsimvec(const arma::mat& z, arma::mat& sigma_log_sim, const arma::vec& variance_intercept,
                  const arma::mat& init, const arma::vec& alpha, const arma::vec& gamma,
                  const arma::vec& beta, const double kappa, const double mu,
                  const arma::ivec& order) {
    const int maxpq = (int) max(order);
    int h = z.n_cols - maxpq; // Assuming z already includes space for burn-in
    int nsim = z.n_rows;
    arma::mat sigma_sim(nsim, h + maxpq);
    arma::mat epsilon(nsim, h + maxpq);
    arma::mat series_sim(nsim, h + maxpq);
    int i,j = 0;
    for(i=maxpq;i<(maxpq+h);i++) {
        sigma_log_sim.col(i).fill(arma::as_scalar(variance_intercept(i)));
        if (order(0) > 0) {
            for(j=0;j<order(0);j++) {
                if((order(0) + j) >= i) {
                    sigma_log_sim.col(i) += alpha(j) * z.col(i - j - 1) + gamma(j) * init.col(j);
                } else {
                    sigma_log_sim.col(i) += alpha(j) * z.col(i - j - 1) + gamma(j) * (abs(z.col(i - j - 1)) - kappa);
                }
            }
        }
        if (order(1) > 0) {
            for(j=0;j<order(1);j++) {
                sigma_log_sim.col(i) += beta(j) * sigma_log_sim.col(i - j - 1);
            }
        }
        sigma_sim.col(i) = sqrt(exp(sigma_log_sim.col(i)));
        epsilon.col(i) = z.col(i) % sigma_sim.col(i);
        series_sim.col(i) = epsilon.col(i) + mu;
    }
    return List::create(Named("sigma") = sigma_sim, Named("series") = series_sim, Named("epsilon") = epsilon);
}

// [[Rcpp::export(.aparchsimvec)]]
List aparchsimvec(arma::mat& epsilon, arma::mat& sigma_power_sim, const arma::mat& z, const arma::vec& variance_intercept,
                  const arma::mat& init, const arma::vec& alpha, const arma::vec& gamma, const arma::vec& beta, const double delta, const double mu,
                  const arma::ivec& order) {
    const int maxpq = (int) max(order);
    int h = z.n_cols - maxpq; // Assuming z already includes space for burn-in
    int nsim = z.n_rows;
    arma::mat series_sim(nsim, h + maxpq);
    arma::mat sigma_sim(nsim, h + maxpq);
    int i,j = 0;
    for(i=maxpq;i<(maxpq+h);i++) {
        sigma_power_sim.col(i).fill(arma::as_scalar(variance_intercept(i)));
        if (order(0) > 0) {
            for(j=0;j<order(0);j++) {
                if((order(0) + j) >= i) {
                    sigma_power_sim.col(i) += alpha(j) * init.col(j);
                } else {
                    sigma_power_sim.col(i) += alpha(j) * pow(abs(epsilon.col(i - j - 1)) - gamma(j) * epsilon.col(i - j - 1), delta);
                }
            }
        }
        if (order(1) > 0) {
            for(j=0;j<order(1);j++) {
                sigma_power_sim.col(i) += beta(j) * sigma_power_sim.col(i - j - 1);
            }
        }
        sigma_sim.col(i) = pow(sigma_power_sim.col(i), 1.0/delta);
        epsilon.col(i) = z.col(i) % sigma_sim.col(i);
        series_sim.col(i) = epsilon.col(i) + mu;
    }
    return List::create(Named("sigma") = sigma_sim, Named("series") = series_sim, Named("epsilon") = epsilon);
}


// [[Rcpp::export(.gjrsimvec)]]
List gjrsimvec(arma::mat& epsilon, arma::mat& sigma_sqr_sim, const arma::mat& z,
               const arma::vec& variance_intercept, const arma::mat& init, const arma::vec& alpha,
               const arma::vec& gamma, const arma::vec& beta, const double mu, const arma::ivec& order)
{
    const int maxpq = (int) max(order);
    int h = z.n_cols - maxpq; // Assuming z already includes space for burn-in
    int nsim = z.n_rows;
    arma::mat sigma_sim(nsim, h + maxpq);
    arma::mat series_sim(nsim, h + maxpq);
    arma::vec tmp(nsim);
    int i,j = 0;
    for(i=maxpq;i<(maxpq+h);i++) {
        sigma_sqr_sim.col(i).fill(arma::as_scalar(variance_intercept(i)));
        if (order(0) > 0) {
            for(j=0;j<order(0);j++) {
                if((order(0) + j) >= i) {
                    sigma_sqr_sim.col(i) += alpha(j) * pow(epsilon.col(i - j - 1), 2.0) + gamma(j) * init.col(j);
                } else {
                    tmp = epsilon.col(i - j - 1);
                    sigma_sqr_sim.col(i) += alpha(j) * pow(epsilon.col(i - j - 1), 2.0) + gamma(j) * (pow(epsilon.col(i - j - 1), 2.0) % (tmp <= 0));
                }
            }
        }
        if (order(1) > 0) {
            for(j=0;j<order(1);j++) {
                sigma_sqr_sim.col(i) += beta(j) * sigma_sqr_sim.col(i - j - 1);
            }
        }
        sigma_sim.col(i) = sqrt(sigma_sqr_sim.col(i));
        epsilon.col(i) = z.col(i) % sigma_sim.col(i);
        series_sim.col(i) = epsilon.col(i) + mu;
    }
    return List::create(Named("sigma") = sigma_sim, Named("series") = series_sim, Named("epsilon") = epsilon);
}

// [[Rcpp::export(.fgarchsimvec)]]
List fgarchsimvec(arma::mat& epsilon, arma::mat& sigma_power_sim, const arma::mat& z, const arma::vec& variance_intercept,
                  const arma::mat& init, const arma::vec& alpha, const arma::vec& gamma, const arma::vec& eta, const arma::vec& beta,
                  const double delta, const double mu, const arma::ivec& order) {
    const int maxpq = (int) max(order);
    int h = z.n_cols - maxpq; // Assuming z already includes space for burn-in
    int nsim = z.n_rows;
    arma::mat series_sim(nsim, h + maxpq);
    arma::mat sigma_sim(nsim, h + maxpq);
    int i,j = 0;
    for(i=maxpq;i<(maxpq+h);i++) {
        sigma_power_sim.col(i).fill(arma::as_scalar(variance_intercept(i)));
        if (order(0) > 0) {
            for(j=0;j<order(0);j++) {
                if((order(0) + j) >= i) {
                    sigma_power_sim.col(i) += alpha(j) * sigma_power_sim.col(i - j - 1) % init.col(j);
                } else {
                    sigma_power_sim.col(i) += alpha(j) * sigma_power_sim.col(i - j - 1) % pow(abs(z.col(i - j - 1) - eta(j)) - gamma(j) * (z.col(i - j - 1) - eta(j)), delta);
                }
            }
        }
        if (order(1) > 0) {
            for(j=0;j<order(1);j++) {
                sigma_power_sim.col(i) += beta(j) * sigma_power_sim.col(i - j - 1);
            }
        }
        sigma_sim.col(i) = pow(sigma_power_sim.col(i), 1.0/delta);
        epsilon.col(i) = z.col(i) % sigma_sim.col(i);
        series_sim.col(i) = epsilon.col(i) + mu;
    }
    return List::create(Named("sigma") = sigma_sim, Named("series") = series_sim, Named("epsilon") = epsilon);
}


// [[Rcpp::export(.cgarchsimvec)]]
List cgarchsimvec(arma::mat& epsilon, arma::mat& sigma_sqr_sim, const arma::mat& z, const arma::vec& variance_intercept,
                  arma::mat& transitory_component_sim, arma::mat& permanent_component_sim, const arma::vec& alpha,
                  const arma::vec& phi, const arma::vec& rho, const arma::vec& beta,
                  const double mu, const arma::ivec& order) {
    const int maxpq = (int) max(order);
    int h = z.n_cols - maxpq; // Assuming z already includes space for burn-in
    int nsim = z.n_rows;
    arma::mat series_sim(nsim, h + maxpq);
    arma::mat sigma_sim(nsim, h + maxpq);
    int i,j = 0;
    for(i=maxpq;i<(maxpq+h);i++) {
        permanent_component_sim.col(i).fill(arma::as_scalar(variance_intercept(i)));
        permanent_component_sim.col(i) += rho(0) * permanent_component_sim.col(i - 1) + phi(0) * (pow(epsilon.col(i - 1), 2.0) - sigma_sqr_sim.col(i - 1));
        if (order(0) > 0) {
            for(j=0;j<order(0);j++) {
                transitory_component_sim.col(i) += alpha(j) * (pow(epsilon.col(i - j - 1), 2.0) - sigma_sqr_sim.col(i - j - 1)) + alpha(j) * transitory_component_sim.col(i - j - 1);
            }
        }
        if (order(1) > 0) {
            for(j=0;j<order(1);j++) {
                transitory_component_sim.col(i) += beta(j) * transitory_component_sim.col(i - j - 1);
            }
        }
        sigma_sqr_sim.col(i) += permanent_component_sim.col(i) + transitory_component_sim.col(i);
        sigma_sim.col(i) = sqrt(sigma_sqr_sim.col(i));
        epsilon.col(i) = z.col(i) % sigma_sim.col(i);
        series_sim.col(i) = epsilon.col(i) + mu;
    }
    return List::create(Named("sigma") = sigma_sim, Named("series") = series_sim, Named("transitory_component") = transitory_component_sim,
                              Named("permanent_component") = permanent_component_sim, Named("epsilon") = epsilon);
}
