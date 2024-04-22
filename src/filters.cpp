#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(.garchfilter)]]
NumericVector garchfilter(const NumericVector residuals, const NumericVector v,
                          const NumericVector initstate, const double omega,
                          const NumericVector alpha, const NumericVector beta,
                          const IntegerVector model) {
    int timesteps = residuals.size();
    // omega [ init_variance, omega]
    // v is the external regressor V x \kappa (already premultiplied in the R code)
    // model [max(p,q) multiplicative ARCH(p) GARCH(q)]
    NumericVector sigma_squared(timesteps);
    NumericVector sigma(timesteps);
    for(int j=0;j<model(0);j++){
        sigma_squared(j) = initstate(j);
    }
    for(int i = model(0);i<timesteps;i++) {
        sigma_squared(i) = omega + v(i);
        if (model(3) > 0.5) sigma_squared(i) = exp(sigma_squared(i));
        if (model(1) > 0) {
            for(int j = 0;j<model(1);j++) {
                sigma_squared(i) += alpha(j) * std::pow(residuals(i - j - 1), 2);
            }
        }
        if (model(2) > 0) {
            for(int j = 0;j<model(2);j++) {
                sigma_squared(i) += beta(j) * sigma_squared(i - j - 1);
            }
        }
    }
    sigma = sqrt(sigma_squared);
    return sigma;
}

// [[Rcpp::export(.egarchfilter)]]
NumericVector egarchfilter(const NumericVector residuals, const NumericVector v,
                           const NumericVector initstate, const double omega,
                           const NumericVector alpha, const NumericVector gamma,
                           const NumericVector beta, const double kappa,
                           const IntegerVector model) {
    int timesteps = residuals.size();
    // omega [ init_variance, omega]
    // v is the external regressor V x \kappa (already premultiplied in the R code)
    // model [max(p,q) multiplicative ARCH(p) GARCH(q)]
    NumericVector sigma_squared(timesteps);
    NumericVector log_sigma_squared(timesteps);
    NumericVector sigma(timesteps);
    NumericVector std_residuals(timesteps);
    for(int j=0;j<model(0);j++){
        sigma_squared(j) = initstate(j);
        log_sigma_squared(j) = log(sigma_squared(j));
        std_residuals(j) = residuals(j)/sqrt(sigma_squared(j));
    }
    for(int i = model(0);i<timesteps;i++) {
        log_sigma_squared(i) = omega + v(i);
        if (model(1) > 0) {
            for(int j = 0;j<model(1);j++) {
                log_sigma_squared(i) += alpha(j) * std_residuals(i - j - 1) + gamma(j) * (fabs(std_residuals(i - j - 1)) - kappa);
            }
        }
        if (model(2) > 0) {
            for(int j = 0;j<model(2);j++) {
                log_sigma_squared(i) += beta(j) * log_sigma_squared(i - j - 1);
            }
        }
        sigma_squared(i) = exp(log_sigma_squared(i));
        sigma(i) = sqrt(sigma_squared(i));
        std_residuals(i) = residuals(i)/sigma(i);
    }
    return sigma;
}


// [[Rcpp::export(.aparchfilter)]]
NumericVector aparchfilter(const NumericVector residuals, const NumericVector v,
                           const NumericVector initstate, const double omega,
                           const NumericVector alpha, const NumericVector gamma,
                           const NumericVector beta,  const double delta,
                           const IntegerVector model) {
    int timesteps = residuals.size();
    // omega [ init_variance, omega]
    // v is the external regressor V x \kappa (already premultiplied in the R code)
    // model [max(p,q) multiplicative ARCH(p) GARCH(q)]
    NumericVector power_sigma(timesteps);
    NumericVector sigma(timesteps);

    for(int j=0;j<model(0);j++){
        sigma(j) = sqrt(initstate(j));
        power_sigma(j) = pow(sigma(j), delta);
    }
    for(int i = model(0);i<timesteps;i++) {
        power_sigma(i) += omega + v(i);
        if (model(3) > 0.5) power_sigma(i) = exp(power_sigma(i));
        if (model(1) > 0) {
            for(int j = 0;j<model(1);j++) {
                power_sigma(i) += alpha(j) * pow(fabs(residuals(i - j - 1)) - gamma(j) * residuals(i - j - 1), delta);
            }
        }
        if (model(2) > 0) {
            for(int j = 0;j<model(2);j++) {
                power_sigma(i) += beta(j) * power_sigma(i - j - 1);
            }
        }
        sigma(i) = pow(power_sigma(i), 1.0/delta);
    }
    return sigma;
}

// [[Rcpp::export(.gjrgarchfilter)]]
NumericVector gjrgarchfilter(const NumericVector residuals, const NumericVector negative_indicator,
                             const NumericVector v, const NumericVector initstate,
                             const double omega, const NumericVector alpha,
                             const NumericVector gamma, const NumericVector beta,
                             const IntegerVector model) {
    int timesteps = residuals.size();
    // omega [ init_variance, omega]
    // v is the external regressor V x \kappa (already premultiplied in the R code)
    // model [max(p,q) multiplicative ARCH(p) GARCH(q)]
    NumericVector sigma_squared(timesteps);
    NumericVector sigma(timesteps);

    for(int j=0;j<model(0);j++){
        sigma(j) = sqrt(initstate(j));
        sigma_squared(j) = initstate(j);
    }
    for(int i = model(0);i<timesteps;i++) {
        sigma_squared(i) += omega + v(i);
        if (model(3) > 0.5) sigma_squared(i) = exp(sigma_squared(i));

        if (model(1) > 0) {
            for(int j = 0;j<model(1);j++) {
                sigma_squared(i) += alpha(j) * pow(residuals(i - j - 1), 2.0) + gamma(j) * (pow(residuals(i - j - 1), 2.0) * negative_indicator(i - j - 1));
            }
        }
        if (model(2) > 0) {
            for(int j = 0;j<model(2);j++) {
                sigma_squared(i) += beta(j) * sigma_squared(i - j - 1);
            }
        }
        sigma(i) = sqrt(sigma_squared(i));
    }
    return sigma;
}

// [[Rcpp::export(.fgarchfilter)]]
NumericVector fgarchfilter(const NumericVector residuals, const NumericVector v,
                           const NumericVector initstate, const double omega,
                           const NumericVector alpha, const NumericVector gamma,
                           const NumericVector eta, const NumericVector beta,
                           const double delta, const IntegerVector model) {
    int timesteps = residuals.size();
    // omega [ init_variance, omega]
    // v is the external regressor V x \kappa (already premultiplied in the R code)
    // model [max(p,q) multiplicative ARCH(p) GARCH(q)]
    NumericVector power_sigma(timesteps);
    NumericVector sigma(timesteps);
    NumericVector std_residuals(timesteps);

    for(int j=0;j<model(0);j++){
        sigma(j) = sqrt(initstate(j));
        power_sigma(j) = pow(sigma(j), delta);
        std_residuals(j) = residuals(j)/sigma(j);
    }
    for(int i = model(0);i<timesteps;i++) {
        power_sigma(i) += omega + v(i);
        if (model(3) > 0.5) power_sigma(i) = exp(power_sigma(i));
        if (model(1) > 0) {
            for(int j = 0;j<model(1);j++) {
                power_sigma(i) += alpha(j) * power_sigma(i - j - 1) * pow(fabs(std_residuals(i - j - 1) - eta(j)) - gamma(j) * (std_residuals(i - j - 1) - eta(j)), delta);
            }
        }
        if (model(2) > 0) {
            for(int j = 0;j<model(2);j++) {
                power_sigma(i) += beta(j) * power_sigma(i - j - 1);
            }
        }
        sigma(i) = pow(power_sigma(i), 1.0/delta);
        std_residuals(i) = residuals(i)/sigma(i);
    }
    return sigma;
}

// [[Rcpp::export(.cgarchfilter)]]
Rcpp::List cgarchfilter(const NumericVector residuals, const NumericVector v,
                        const NumericMatrix initstate, const double omega,
                        const NumericVector alpha, const NumericVector rho,
                        const NumericVector phi, const NumericVector beta,
                        const IntegerVector model) {
    int timesteps = residuals.size();
    // omega [ omega]
    // v is the external regressor V x \kappa (already premultiplied in the R code)
    // model [max(p,q) multiplicative ARCH(p) GARCH(q)]
    NumericVector sigma(timesteps);
    NumericVector sigma_squared(timesteps);
    NumericVector permanent_component(timesteps);
    NumericVector transitory_component(timesteps);
    NumericVector residuals_squared = pow(residuals, 2.0);
    int j,i;
    for(j=0;j<model(0);j++){
        double sinit = initstate(j,0);
        double pinit = initstate(j,1);
        permanent_component(j) = pinit;
        transitory_component(j) = sinit;
        sigma_squared(j) = pinit + sinit;
        sigma(j) = sqrt(sigma_squared(j));
    }
    for(i = model(0);i<timesteps;i++) {
        permanent_component(i) += omega + v(i);
        if (model(3) > 0.5) permanent_component(i) = exp(permanent_component(i));
        permanent_component(i) += rho(0) * permanent_component(i - 1) + phi(0) * (residuals_squared(i - 1) - sigma_squared(i - 1));
        if (model(1) > 0) {
            for(j = 0;j<model(1);j++) {
                transitory_component(i) += alpha(j) * transitory_component(i - j - 1) + alpha(j) *  (residuals_squared(i - j - 1) - sigma_squared(i - j - 1));
            }
        }
        if (model(2) > 0) {
            for(j = 0;j<model(2);j++) {
                transitory_component(i) += beta(j) *  transitory_component(i - j - 1);
            }
        }
        sigma_squared(i) = permanent_component(i) + transitory_component(i);
        sigma(i) += sqrt(sigma_squared(i));
    }
    Rcpp::List output = Rcpp::List::create(Rcpp::Named("sigma") = wrap(sigma),
                                           Rcpp::Named("transitory_component") = wrap(transitory_component),
                                           Rcpp::Named("permanent_component") = wrap(permanent_component));
    return output;
}
