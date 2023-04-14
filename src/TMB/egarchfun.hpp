/// @file egarchfun.hpp
#ifndef egarchfun_hpp
#define egarchfun_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type egarchfun(objective_function<Type>* obj) {
    DATA_VECTOR(y);
    // variance initialization
    DATA_SCALAR(backcast_lambda);
    DATA_INTEGER(samplen);
    DATA_STRING(initmethod);
    PARAMETER(mu);
    PARAMETER(omega);
    PARAMETER_VECTOR(alpha);
    PARAMETER_VECTOR(gamma);
    PARAMETER_VECTOR(beta);
    PARAMETER_VECTOR(xi);
    PARAMETER_VECTOR(distribution);
    // parameter scaling vector
    DATA_VECTOR(pscale);
    // variance regressors
    DATA_MATRIX(v);
    // model flags
    DATA_IVECTOR(cmodel);
    const int timesteps = y.rows();
    vector<Type> regressors(timesteps);
    vector<Type> std_residuals(timesteps);
    vector<Type> sigma_squared(timesteps);
    vector<Type> sigma(timesteps);
    vector<Type> log_sigma_squared(timesteps);
    sigma_squared.setZero();
    sigma.setZero();
    std_residuals.setZero();
    regressors.setZero();
    log_sigma_squared.setZero();
    int dclass = cmodel(5);
    int m = v.cols();
    int j = 0;
    // re-scale parameters
    int k = 0;
    mu *= pscale(k);
    k += 1;
    omega *= pscale(k);
    k += 1;
    for(j = 0;j<cmodel(1);j++) {
        alpha(j) *= pscale(j + k);
    }
    if (cmodel(1) == 0){
        k += 1;
    } else {
        k += cmodel(1);
    }
    for(j = 0;j<cmodel(1);j++){
        gamma(j) *= pscale(j + k);
    }
    if (cmodel(1) == 0){
        k += 1;
    } else {
        k += cmodel(1);
    }
    for(j = 0;j<cmodel(2);j++){
        beta(j) *= pscale(j + k);
    }
    if (cmodel(2) == 0){
        k += 1;
    } else {
        k += cmodel(2);
    }
    for(j = 0;j<m;j++){
        xi(j) *= pscale(j + k);
    }
    k += m;
    distribution(0) *= pscale(k);
    distribution(1) *= pscale(k + 1);
    distribution(2) *= pscale(k + 2);
    vector<Type> residuals = y.array() - mu;
    vector<Type> tmp_block = residuals.tail(timesteps - cmodel(0));
    Type initial_variance = garchextra::init_power_variance(tmp_block, initmethod, backcast_lambda, Type(2.0), samplen);
    Type initial_log_variance = log(initial_variance);
    //Type initial_variance = exp(initial_log_variance);
    Type skew = distribution(0);
    Type shape = distribution(1);
    Type lambda = distribution(2);
    Type kappa = egarchkappa::egarch_moment_func(skew, shape, lambda, dclass);
    ADREPORT(kappa);
    // initialize variance
    for(j = 0;j<cmodel(0);j++) {
        sigma_squared(j) += initial_variance;
        log_sigma_squared(j) = initial_log_variance;
        // zero out the initial values
        residuals(j) = 0.0;
        sigma(j) = sqrt(sigma_squared(j));
    }
    Type persistence = beta.sum();
    regressors = v * xi;
    vector<Type> residuals_squared = residuals.array().square();
    vector<Type> variance_intercept(timesteps);
    Type target_omega = 0.0;
    if (cmodel(3) > 0.5) {
        Type log_sample_variance = log(residuals.tail(timesteps - cmodel(0)).square().mean());
        target_omega = log_sample_variance * (Type(1.0) - persistence);
        // subtract (mean of v) * xi
        vector<Type> meanc = v.colwise().mean();
        Type mean_regressors = (meanc.array() * xi.array()).sum();
        target_omega -= mean_regressors;
        variance_intercept.fill(target_omega);
        ADREPORT(target_omega);
    } else {
        target_omega = omega;
        variance_intercept.fill(omega);
        ADREPORT(target_omega);
    }
    ADREPORT(persistence);
    variance_intercept.array() = variance_intercept.array() + regressors.array();

    for(int i = cmodel(0);i<timesteps;i++){
        log_sigma_squared(i) += variance_intercept(i);
        for(j = 0;j<cmodel(1);j++){
            log_sigma_squared(i) += alpha(j) * std_residuals(i - j - 1) + gamma(j) * (fabs(std_residuals(i - j - 1)) - kappa);
        }
        for(j = 0;j<cmodel(2);j++){
            log_sigma_squared(i) += beta(j) * log_sigma_squared(i - j - 1);
        }
        sigma_squared(i) = exp(log_sigma_squared(i));
        sigma(i) = sqrt(sigma_squared(i));
        std_residuals(i) = residuals(i)/sigma(i);
    }
    vector<Type> tmp_vector = distfun::distlike(std_residuals, distribution(0), distribution(1), distribution(2), dclass)/sigma.array();
    vector<Type> ll_vector = tmp_vector.tail(timesteps - cmodel(0));
    REPORT(target_omega);
    REPORT(kappa);
    REPORT(initial_variance);
    REPORT(sigma);
    REPORT(ll_vector);
    Type nll = Type(-1.0) * ll_vector.log().sum();
    return(nll);
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
