#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;

// [[Rcpp::export(.garchsimvec)]]
List garchsimvec(Eigen::Map<Eigen::MatrixXd>& epsilon, Eigen::Map<Eigen::MatrixXd>& sigma_sqr_sim, const Eigen::Map<Eigen::MatrixXd>& z,
                 const Eigen::Map<Eigen::VectorXd>& variance_intercept, const Eigen::Map<Eigen::MatrixXd>& init, const Eigen::Map<Eigen::VectorXd>& alpha,
                 const Eigen::Map<Eigen::VectorXd>& beta, const double mu, const Eigen::Map<Eigen::VectorXi>& order) {
    const int maxpq = (int) order.maxCoeff();
    int h = z.cols() - maxpq; // Assuming z already includes space for burn-in
    int nsim = z.rows();
    Eigen::MatrixXd sigma_sim = Eigen::MatrixXd::Zero(nsim, h + maxpq);
    Eigen::MatrixXd series_sim = Eigen::MatrixXd::Zero(nsim, h + maxpq);
    int i,j = 0;
    for(i=maxpq; i<(maxpq+h); i++) {
        sigma_sqr_sim.col(i).setConstant(variance_intercept(i));
        if (order(0) > 0) {
            for(j=0; j<order(0); j++) {
                if((order(0) + j) >= i) {
                    sigma_sqr_sim.col(i) += alpha(j) * init.col(j);
                } else {
                    sigma_sqr_sim.col(i) += alpha(j) * (epsilon.col(i - j - 1).array().square()).matrix();
                }
            }
        }
        if (order(1) > 0) {
            for(j=0; j<order(1); j++) {
                sigma_sqr_sim.col(i) += beta(j) * sigma_sqr_sim.col(i - j - 1);
            }
        }
        sigma_sim.col(i) = sigma_sqr_sim.col(i).cwiseSqrt();
        epsilon.col(i) = z.col(i).array() * sigma_sim.col(i).array();
        series_sim.col(i) = epsilon.col(i).array() + mu;
    }    return List::create(Named("sigma") = sigma_sim, Named("series") = series_sim, Named("epsilon") = epsilon);
}

// [[Rcpp::export(.egarchsimvec)]]
List egarchsimvec(const Eigen::MatrixXd& z, Eigen::MatrixXd& sigma_log_sim, const Eigen::VectorXd& variance_intercept,
                  const Eigen::MatrixXd& init, const Eigen::VectorXd& alpha, const Eigen::VectorXd& gamma,
                  const Eigen::VectorXd& beta, const double kappa, const double mu,
                  const Eigen::VectorXi& order) {
    const int maxpq = order.maxCoeff();
    int h = z.cols() - maxpq;
    int nsim = z.rows();
    Eigen::MatrixXd sigma_sim = Eigen::MatrixXd::Zero(nsim, h + maxpq);
    Eigen::MatrixXd epsilon = Eigen::MatrixXd::Zero(nsim, h + maxpq);
    Eigen::MatrixXd series_sim = Eigen::MatrixXd::Zero(nsim, h + maxpq);
    int i, j;
    for(i = maxpq; i < (maxpq + h); i++) {
        sigma_log_sim.col(i).setConstant(variance_intercept(i));
        if (order(0) > 0) {
            for(j = 0; j < order(0); j++) {
                if((order(0) + j) >= i) {
                    sigma_log_sim.col(i) += alpha(j) * z.col(i - j - 1) + gamma(j) * init.col(j);
                } else {
                    sigma_log_sim.col(i).array() += alpha(j) * z.col(i - j - 1).array() + gamma(j) * (z.col(i - j - 1).array().abs() - kappa);
                }
            }
        }
        if (order(1) > 0) {
            for(j = 0; j < order(1); j++) {
                sigma_log_sim.col(i) += beta(j) * sigma_log_sim.col(i - j - 1);
            }
        }
        sigma_sim.col(i) = (sigma_log_sim.col(i).array().exp()).sqrt();
        epsilon.col(i) = z.col(i).array() * sigma_sim.col(i).array();
        series_sim.col(i) = epsilon.col(i).array() + mu;
    }
    return List::create(Named("sigma") = sigma_sim, Named("series") = series_sim, Named("epsilon") = epsilon);
}

// [[Rcpp::export(.aparchsimvec)]]
List aparchsimvec(Eigen::MatrixXd& epsilon, Eigen::MatrixXd& sigma_power_sim, const Eigen::MatrixXd& z,
                   const Eigen::VectorXd& variance_intercept, const Eigen::MatrixXd& init,
                   const Eigen::VectorXd& alpha, const Eigen::VectorXd& gamma, const Eigen::VectorXd& beta,
                   const double delta, const double mu, const Eigen::VectorXi& order) {
    const int maxpq = order.maxCoeff();
    int h = z.cols() - maxpq; // Assuming z already includes space for burn-in
    int nsim = z.rows();
    Eigen::MatrixXd series_sim = Eigen::MatrixXd::Zero(nsim, h + maxpq);
    Eigen::MatrixXd sigma_sim = Eigen::MatrixXd::Zero(nsim, h + maxpq);
    int i, j;
    for (i = maxpq; i < (maxpq + h); i++) {
        sigma_power_sim.col(i).setConstant(variance_intercept(i));
        if (order(0) > 0) {
            for (j = 0; j < order(0); j++) {
                if ((order(0) + j) >= i) {
                    sigma_power_sim.col(i) += alpha(j) * init.col(j);
                } else {
                    sigma_power_sim.col(i).array() += alpha(j) *
                        pow((epsilon.col(i - j - 1).array().abs() - gamma(j) * epsilon.col(i - j - 1).array()), delta);
                }
            }
        }
        if (order(1) > 0) {
            for (j = 0; j < order(1); j++) {
                sigma_power_sim.col(i) += beta(j) * sigma_power_sim.col(i - j - 1);
            }
        }
        sigma_sim.col(i) = sigma_power_sim.col(i).array().pow(1.0 / delta);
        epsilon.col(i) = z.col(i).array() * sigma_sim.col(i).array();
        series_sim.col(i) = epsilon.col(i).array() + mu;
    }
    return List::create(Named("sigma") = sigma_sim, Named("series") = series_sim, Named("epsilon") = epsilon);
}

// [[Rcpp::export(.gjrsimvec)]]
List gjrsimvec(Eigen::MatrixXd& epsilon, Eigen::MatrixXd& sigma_sqr_sim, const Eigen::MatrixXd& z,
                const Eigen::VectorXd& variance_intercept, const Eigen::MatrixXd& init, const Eigen::VectorXd& alpha,
                const Eigen::VectorXd& gamma, const Eigen::VectorXd& beta, const double mu, const Eigen::VectorXi& order)
{
    const int maxpq = order.maxCoeff();
    int h = z.cols() - maxpq; // Assuming z already includes space for burn-in
    int nsim = z.rows();
    Eigen::MatrixXd sigma_sim = Eigen::MatrixXd::Zero(nsim, h + maxpq);
    Eigen::MatrixXd series_sim = Eigen::MatrixXd::Zero(nsim, h + maxpq);
    Eigen::VectorXd tmp = Eigen::VectorXd::Zero(nsim);
    int i, j = 0;
    for(i = maxpq; i < (maxpq + h); i++) {
        sigma_sqr_sim.col(i).setConstant(variance_intercept(i));
        if (order(0) > 0) {
            for(j = 0; j < order(0); j++) {
                if((order(0) + j) >= i) {
                    sigma_sqr_sim.col(i).array() += alpha(j) * epsilon.col(i - j - 1).array().square() + gamma(j) * init.col(j).array();
                } else {
                    tmp = epsilon.col(i - j - 1);
                    sigma_sqr_sim.col(i).array() += alpha(j) * epsilon.col(i - j - 1).array().square() +
                        gamma(j) * (epsilon.col(i - j - 1).array().square() * (tmp.array() <= 0).cast<double>());
                }
            }
        }
        if (order(1) > 0) {
            for(j = 0; j < order(1); j++) {
                sigma_sqr_sim.col(i) += beta(j) * sigma_sqr_sim.col(i - j - 1);
            }
        }
        sigma_sim.col(i) = sigma_sqr_sim.col(i).array().sqrt();
        epsilon.col(i) = z.col(i).array() * sigma_sim.col(i).array();
        series_sim.col(i) = epsilon.col(i).array() + mu;
    }
    return List::create(Named("sigma") = sigma_sim, Named("series") = series_sim, Named("epsilon") = epsilon);
}

// [[Rcpp::export(.fgarchsimvec)]]
List fgarchsimvec(Eigen::MatrixXd& epsilon, Eigen::MatrixXd& sigma_power_sim, const Eigen::MatrixXd& z,
                   const Eigen::VectorXd& variance_intercept, const Eigen::MatrixXd& init, const Eigen::VectorXd& alpha,
                   const Eigen::VectorXd& gamma, const Eigen::VectorXd& eta, const Eigen::VectorXd& beta,
                   const double delta, const double mu, const Eigen::VectorXi& order) {
    const int maxpq = order.maxCoeff();
    int h = z.cols() - maxpq; // Assuming z already includes space for burn-in
    int nsim = z.rows();
    Eigen::MatrixXd series_sim = Eigen::MatrixXd::Zero(nsim, h + maxpq);
    Eigen::MatrixXd sigma_sim = Eigen::MatrixXd::Zero(nsim, h + maxpq);
    int i, j;
    for (i = maxpq; i < (maxpq + h); i++) {
        sigma_power_sim.col(i).setConstant(variance_intercept(i));
        if (order(0) > 0) {
            for (j = 0; j < order(0); j++) {
                if ((order(0) + j) >= i) {
                    sigma_power_sim.col(i).array() += alpha(j) * sigma_power_sim.col(i - j - 1).array() * init.col(j).array();
                } else {
                    sigma_power_sim.col(i).array() += alpha(j) * sigma_power_sim.col(i - j - 1).array() *
                        (abs((z.col(i - j - 1).array() - eta(j)).array()) - gamma(j) * (z.col(i - j - 1).array() - eta(j)).array()).pow(delta);
                }
            }
        }
        if (order(1) > 0) {
            for (j = 0; j < order(1); j++) {
                sigma_power_sim.col(i) += beta(j) * sigma_power_sim.col(i - j - 1);
            }
        }
        sigma_sim.col(i) = sigma_power_sim.col(i).array().pow(1.0 / delta);
        epsilon.col(i) = z.col(i).array() * sigma_sim.col(i).array();
        series_sim.col(i) = epsilon.col(i).array() + mu;
    }
    return List::create(Named("sigma") = sigma_sim, Named("series") = series_sim, Named("epsilon") = epsilon);
}

// [[Rcpp::export(.cgarchsimvec)]]
List cgarchsimvec(Eigen::MatrixXd& epsilon, Eigen::MatrixXd& sigma_sqr_sim, const Eigen::MatrixXd& z,
                   const Eigen::VectorXd& variance_intercept, Eigen::MatrixXd& transitory_component_sim,
                   Eigen::MatrixXd& permanent_component_sim, const Eigen::VectorXd& alpha,
                   const Eigen::VectorXd& phi, const Eigen::VectorXd& rho, const Eigen::VectorXd& beta,
                   const double mu, const Eigen::VectorXi& order) {
    const int maxpq = order.maxCoeff();
    int h = z.cols() - maxpq; // Assuming z already includes space for burn-in
    int nsim = z.rows();
    Eigen::MatrixXd series_sim = Eigen::MatrixXd::Zero(nsim, h + maxpq);
    Eigen::MatrixXd sigma_sim = Eigen::MatrixXd::Zero(nsim, h + maxpq);
    int i, j;
    for(i = maxpq; i < (maxpq + h); i++) {
        permanent_component_sim.col(i).setConstant(variance_intercept(i));
        if (i > 0) { // Ensure i-1 is valid
            permanent_component_sim.col(i).array() += rho(0) * permanent_component_sim.col(i - 1).array() +
                phi(0) * (epsilon.col(i - 1).array().square() - sigma_sqr_sim.col(i - 1).array());
        }
        if (order(0) > 0) {
            for(j = 0; j < order(0); j++) {
                if(i - j - 1 >= 0) { // Ensure the index is valid
                    transitory_component_sim.col(i).array() += alpha(j) *
                        (epsilon.col(i - j - 1).array().square() - sigma_sqr_sim.col(i - j - 1).array()) +
                        alpha(j) * transitory_component_sim.col(i - j - 1).array();
                }
            }
        }
        if (order(1) > 0) {
            for(j = 0; j < order(1); j++) {
                if(i - j - 1 >= 0) { // Ensure the index is valid
                    transitory_component_sim.col(i) += beta(j) * transitory_component_sim.col(i - j - 1);
                }
            }
        }
        sigma_sqr_sim.col(i) += permanent_component_sim.col(i) + transitory_component_sim.col(i);
        sigma_sim.col(i) = sigma_sqr_sim.col(i).array().sqrt();
        epsilon.col(i) = z.col(i).array() * sigma_sim.col(i).array();
        series_sim.col(i) = epsilon.col(i).array() + mu;
    }
    return List::create(Named("sigma") = sigma_sim, Named("series") = series_sim, Named("transitory_component") = transitory_component_sim,
                              Named("permanent_component") = permanent_component_sim, Named("epsilon") = epsilon);
}
