#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(.garchfilter)]]
NumericVector garchfilter(NumericVector residuals, NumericVector v, NumericVector initstate, NumericVector omega, NumericVector alpha, NumericVector beta, IntegerVector model) {
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
        sigma_squared(i) = omega(0) + v(i);
        if (model(1) > 0.5) sigma_squared(i) = exp(sigma_squared(i));
        if (model(2) > 0) {
            for(int j = 0;j<model(2);j++) {
                sigma_squared(i) += alpha(j) * std::pow(residuals(i - j - 1), 2);
            }
        }
        if (model(3) > 0) {
            for(int j = 0;j<model(3);j++) {
                sigma_squared(i) += beta(j) * sigma_squared(i - j - 1);
            }
        }
    }
    sigma = sqrt(sigma_squared);
    return sigma;
}

// [[Rcpp::export(.egarchfilter)]]
NumericVector egarchfilter(NumericVector residuals, NumericVector v, NumericVector initstate,
                           NumericVector omega, NumericVector alpha, NumericVector gamma,
                           NumericVector beta, NumericVector kappa, IntegerVector model) {
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
        std_residuals(j) = 0.0;
    }
    for(int i = model(0);i<timesteps;i++) {
        log_sigma_squared(i) = omega(0) + v(i);
        if (model(2) > 0) {
            for(int j = 0;j<model(2);j++) {
                log_sigma_squared(i) += alpha(j) * std_residuals(i - j - 1) + gamma(j) * (abs(std_residuals(i - j - 1)) - kappa(0));
            }
        }
        if (model(3) > 0) {
            for(int j = 0;j<model(3);j++) {
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
NumericVector aparchfilter(NumericVector residuals, NumericVector v, NumericVector initstate,
                           NumericVector omega, NumericVector alpha, NumericVector gamma,
                           NumericVector beta,  NumericVector delta, IntegerVector model) {
    int timesteps = residuals.size();
    // omega [ init_variance, omega]
    // v is the external regressor V x \kappa (already premultiplied in the R code)
    // model [max(p,q) multiplicative ARCH(p) GARCH(q)]
    NumericVector power_sigma(timesteps);
    NumericVector sigma(timesteps);

    for(int j=0;j<model(0);j++){
        sigma(j) = sqrt(initstate(j));
        power_sigma(j) = pow(sigma(j), delta(0));
    }
    for(int i = model(0);i<timesteps;i++) {
        power_sigma(i) += omega(0) + v(i);
        if (model(2) > 0) {
            for(int j = 0;j<model(2);j++) {
                power_sigma(i) += alpha(j) * pow(fabs(residuals(i - j - 1)) - gamma(j) * residuals(i - j - 1), delta(0));
            }
        }
        if (model(3) > 0) {
            for(int j = 0;j<model(3);j++) {
                power_sigma(i) += beta(j) * power_sigma(i - j - 1);
            }
        }
        sigma(i) = pow(power_sigma(i), 1.0/delta(0));
    }
    return sigma;
}

// [[Rcpp::export(.gjrgarchfilter)]]
NumericVector gjrgarchfilter(NumericVector residuals, NumericVector negative_indicator, NumericVector v,
                             NumericVector initstate, NumericVector omega, NumericVector alpha, NumericVector gamma,
                             NumericVector beta, IntegerVector model) {
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
        sigma_squared(i) += omega(0) + v(i);
        if (model(2) > 0) {
            for(int j = 0;j<model(2);j++) {
                sigma_squared(i) += alpha(j) * pow(residuals(i - j - 1), 2.0) + gamma(j) * (pow(residuals(i - j - 1), 2.0) * negative_indicator(i - j - 1));
            }
        }
        if (model(3) > 0) {
            for(int j = 0;j<model(3);j++) {
                sigma_squared(i) += beta(j) * sigma_squared(i - j - 1);
            }
        }
        sigma(i) = sqrt(sigma_squared(i));
    }
    return sigma;
}

// [[Rcpp::export(.fgarchfilter)]]
NumericVector fgarchfilter(NumericVector residuals, NumericVector v, NumericVector initstate,
                           NumericVector omega, NumericVector alpha, NumericVector gamma,
                           NumericVector eta, NumericVector beta,  NumericVector delta,
                           IntegerVector model) {
    int timesteps = residuals.size();
    // omega [ init_variance, omega]
    // v is the external regressor V x \kappa (already premultiplied in the R code)
    // model [max(p,q) multiplicative ARCH(p) GARCH(q)]
    NumericVector power_sigma(timesteps);
    NumericVector sigma(timesteps);
    NumericVector std_residuals(timesteps);

    for(int j=0;j<model(0);j++){
        sigma(j) = sqrt(initstate(j));
        power_sigma(j) = pow(sigma(j), delta(0));
        std_residuals(j) = 0.0;
    }
    for(int i = model(0);i<timesteps;i++) {
        power_sigma(i) += omega(0) + v(i);
        if (model(2) > 0) {
            for(int j = 0;j<model(2);j++) {
                power_sigma(i) += alpha(j) * power_sigma(i - j - 1) * pow(fabs(std_residuals(i - j - 1) - eta(j)) - gamma(j) * (std_residuals(i - j - 1) - eta(j)), delta(0));
            }
        }
        if (model(3) > 0) {
            for(int j = 0;j<model(3);j++) {
                power_sigma(i) += beta(j) * power_sigma(i - j - 1);
            }
        }
        sigma(i) = pow(power_sigma(i), 1.0/delta(0));
        std_residuals(i) = residuals(i)/sigma(i);
    }
    return sigma;
}

// [[Rcpp::export(.cgarchfilter)]]
Rcpp::List cgarchfilter(NumericVector residuals, NumericVector v, NumericMatrix initstate,
                           NumericVector omega, NumericVector alpha, NumericVector rho,
                           NumericVector phi, NumericVector beta, IntegerVector model) {
    int timesteps = residuals.size();
    // omega [ omega]
    // v is the external regressor V x \kappa (already premultiplied in the R code)
    // model [max(p,q) multiplicative ARCH(p) GARCH(q)]
    NumericVector sigma(timesteps);
    NumericVector sigma_squared(timesteps);
    NumericVector permanent_component(timesteps);
    NumericVector residuals_squared = pow(residuals, 2.0);
    int j,i;
    for(j=0;j<model(0);j++){
        double pinit = initstate(j,1);
        double vinit = initstate(j,0);
        if (model(1) > 0.5) {
            pinit = exp(pinit);
        }
        permanent_component(j) = pinit;
        sigma_squared(j) += vinit + permanent_component(j);
        sigma(j) += sqrt(sigma_squared(j));
        residuals_squared(j) = 0.0;
    }
    for(i = model(0);i<timesteps;i++) {
        permanent_component(i) += omega(0) + v(i);
        if (model(1) > 0.5) permanent_component(i) = exp(permanent_component(i));
        permanent_component(i) += rho(0) * permanent_component(i - 1) + phi(0) * (residuals_squared(i - 1) - sigma_squared(i - 1));
        sigma_squared(i) += permanent_component(i);
        if (model(2) > 0) {
            for(j = 0;j<model(2);j++) {
                sigma_squared(i) += alpha(j) * (residuals_squared(i - j - 1) - permanent_component(i - j - 1));
            }
        }
        if (model(3) > 0) {
            for(j = 0;j<model(3);j++) {
                sigma_squared(i) += beta(j) * (sigma_squared(i - j - 1) - permanent_component(i - j - 1));
            }
        }
        sigma(i) += sqrt(sigma_squared(i));
    }
    Rcpp::List output = Rcpp::List::create(Rcpp::Named("sigma") = wrap(sigma),
                                           Rcpp::Named("permanent_component") = wrap(permanent_component));
    return output;
}
