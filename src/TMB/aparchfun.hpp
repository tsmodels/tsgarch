/// @file aparchfun.hpp
#ifndef aparchfun_hpp
#define aparchfun_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type aparchfun(objective_function<Type>* obj) {
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
    PARAMETER(delta);
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
    vector<Type> sigma_power(timesteps);
    vector<Type> sigma(timesteps);
    vector<Type> std_residuals(timesteps);
    sigma_power.setZero();
    sigma.setZero();
    std_residuals.setZero();
    regressors.setZero();
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
    delta *= pscale(k);
    k += 1;
    for(j = 0;j<m;j++){
        xi(j) *= pscale(j + k);
    }
    k += m;
    distribution(0) *= pscale(k);
    distribution(1) *= pscale(k + 1);
    distribution(2) *= pscale(k + 2);
    vector<Type> residuals = y.array() - mu;
    vector<Type> tmp_block = residuals.tail(timesteps - cmodel(0));
    Type initial_power_sigma = garchextra::init_power_variance(tmp_block, initmethod, backcast_lambda, delta, samplen);
    Type initial_variance = pow(initial_power_sigma, Type(2.0)/delta);
    // initialize variance
    for(j = 0;j<cmodel(0);j++) {
        sigma_power(j) += initial_power_sigma;
        // zero out the initial values
        residuals(j) = 0.0;
        sigma(j) = pow(sigma_power(j), Type(1.0)/delta);
    }
    Type persistence = beta.sum();
    vector<Type> kappa(cmodel(1));
    for(j = 0;j<cmodel(1);j++) {
        kappa(j) = aparchkappa::aparch_moment_func(gamma(j), delta, distribution(0), distribution(1), distribution(2), dclass);
        persistence += alpha(j) * kappa(j);
    }
    ADREPORT(persistence);
    ADREPORT(kappa);
    regressors = v * xi;
    vector<Type> power_sigma_intercept(timesteps);
    Type target_omega = 0.0;
    if (cmodel(3) > 0.5) {
        Type sample_power_sigma = residuals.tail(timesteps - cmodel(0)).abs().pow(delta).mean();
        target_omega = sample_power_sigma * (Type(1.0) - persistence);
        // subtract (mean of v) * xi
        vector<Type> meanc = v.colwise().mean();
        Type mean_regressors = (meanc.array() * xi.array()).sum();
        target_omega -= mean_regressors;
        power_sigma_intercept.fill(target_omega);
        ADREPORT(target_omega);
    } else {
        target_omega = omega;
        power_sigma_intercept.fill(omega);
        ADREPORT(target_omega);
    }
    power_sigma_intercept.array() = power_sigma_intercept.array() + regressors.array();
    if (cmodel(4) > 0.5) power_sigma_intercept = power_sigma_intercept.array().exp();

    for(int i = cmodel(0);i<timesteps;i++){
        sigma_power(i) += power_sigma_intercept(i);
        for(j = 0;j<cmodel(1);j++){
            sigma_power(i) += alpha(j) * pow(fabs(residuals(i - j - 1)) - gamma(j) * residuals(i - j - 1), delta);
        }
        for(j = 0;j<cmodel(2);j++){
            sigma_power(i) += beta(j) * sigma_power(i - j - 1);
        }
        sigma(i) = pow(sigma_power(i), Type(1.0)/delta);
        std_residuals(i) = residuals(i)/sigma(i);
    }
    vector<Type> tmp_vector = distfun::distlike(std_residuals, distribution(0), distribution(1), distribution(2), dclass)/sigma.array();
    vector<Type> ll_vector = tmp_vector.tail(timesteps - cmodel(0));
    REPORT(target_omega);
    REPORT(persistence);
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
