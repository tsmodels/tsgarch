#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;

// [[Rcpp::export(.garchsimvec)]]
List garchsimvec(Eigen::Map<Eigen::MatrixXd>& epsilon, Eigen::Map<Eigen::MatrixXd>& sigma_sqr_sim, const Eigen::Map<Eigen::MatrixXd>& z,
                 const Eigen::Map<Eigen::VectorXd>& variance_intercept, const Eigen::Map<Eigen::MatrixXd>& init, const Eigen::Map<Eigen::VectorXd>& alpha,
                 const Eigen::Map<Eigen::VectorXd>& beta, const double mu, const Eigen::Map<Eigen::VectorXi>& order, const int presample) {
    // presample is the pre-sample/burn-in column count to use, which the R
    // wrapper sets to max(garch order, arma order) so that this recursion
    // and the (optional) ARMA mean overlay in .armasimvec() below share the
    // same presample boundary (mirroring rugarch's combined maxOrder used
    // by both its variance (sgarchsimC) and mean (armaxsim) simulation
    // routines). The ARCH/GARCH lag lookback logic below is unaffected by
    // how large the presample region is, since it always indexes backward
    // by an absolute column offset from i.
    const int maxpq = presample;
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

// ARMA mean equation overlay for simulated paths, applied on top of an
// already-simulated GARCH variance/innovation process (epsilon = sigma * z,
// from .garchsimvec() above, unmodified). This mirrors rugarch's own
// two-stage design (variance simulated by sgarchsimC, then the ARMA mean
// overlaid by armaxsim using those same innovations - see
// rugarch's src/garchsim.cpp: msgarchsim() + marmaxsim()):
//   x_t = mu + sum_i ar_i * (x_{t-i} - mu) + sum_j ma_j * eps_{t-j} + eps_t
// AR feedback uses the (already-simulated, mu-seeded) series itself so that
// the joint stochastic dependence between the mean and variance processes
// is preserved; MA feedback and the current-period shock use the same
// eps_t driving the variance recursion, exactly as in the joint ARMA-GARCH
// data-generating process. series_sim's pre-sample columns (< presample)
// must be seeded to mu by the caller so that the AR lookback term vanishes
// there (x_{t-i} - mu = 0), matching a "start from the unconditional mean"
// convention.
// [[Rcpp::export(.armasimvec)]]
Eigen::MatrixXd armasimvec(Eigen::Map<Eigen::MatrixXd>& series_sim, const Eigen::Map<Eigen::MatrixXd>& epsilon,
                           const Eigen::Map<Eigen::VectorXd>& ar, const Eigen::Map<Eigen::VectorXd>& ma,
                           const double mu, const int presample) {
    const int T = (int) series_sim.cols();
    const int ar_order = (int) ar.size();
    const int ma_order = (int) ma.size();
    int i, j;
    for (i = presample; i < T; i++) {
        series_sim.col(i).setConstant(mu);
        for (j = 0; j < ar_order; j++) {
            series_sim.col(i) += ar(j) * (series_sim.col(i - j - 1).array() - mu).matrix();
        }
        for (j = 0; j < ma_order; j++) {
            series_sim.col(i) += ma(j) * epsilon.col(i - j - 1);
        }
        series_sim.col(i) += epsilon.col(i);
    }
    return series_sim;
}

// [[Rcpp::export(.egarchsimvec)]]
List egarchsimvec(const Eigen::MatrixXd& z, Eigen::MatrixXd& sigma_log_sim, const Eigen::VectorXd& variance_intercept,
                  const Eigen::MatrixXd& init, const Eigen::VectorXd& alpha, const Eigen::VectorXd& gamma,
                  const Eigen::VectorXd& beta, const double kappa, const double mu,
                  const Eigen::VectorXi& order, const int presample) {
    // see .garchsimvec() above for why this is an explicit parameter rather
    // than derived from order.maxCoeff()
    const int maxpq = presample;
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
                   const double delta, const double mu, const Eigen::VectorXi& order, const int presample) {
    // see .garchsimvec() above for why this is an explicit parameter rather
    // than derived from order.maxCoeff()
    const int maxpq = presample;
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
                const Eigen::VectorXd& gamma, const Eigen::VectorXd& beta, const double mu, const Eigen::VectorXi& order, const int presample)
{
    // see .garchsimvec() above for why this is an explicit parameter rather
    // than derived from order.maxCoeff()
    const int maxpq = presample;
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
                   const double delta, const double mu, const Eigen::VectorXi& order, const int presample) {
    // see .garchsimvec() above for why this is an explicit parameter rather
    // than derived from order.maxCoeff()
    const int maxpq = presample;
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
                   const double mu, const Eigen::VectorXi& order, const int presample) {
    // see .garchsimvec() above for why this is an explicit parameter rather
    // than derived from order.maxCoeff()
    const int maxpq = presample;
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
