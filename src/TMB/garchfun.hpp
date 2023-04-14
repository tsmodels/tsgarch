/// @file garchfun.hpp
#ifndef garchfun_hpp
#define garchfun_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type garchfun(objective_function<Type>* obj) {
    DATA_VECTOR(y);
    // variance initialization
    DATA_SCALAR(backcast_lambda);
    DATA_INTEGER(samplen);
    DATA_STRING(initmethod);
    PARAMETER(mu);
    PARAMETER(omega);
    PARAMETER_VECTOR(alpha);
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
    vector<Type> sigma_squared(timesteps);
    sigma_squared.setZero();
    regressors.setZero();
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
    // initialize variance
    for(j = 0;j<cmodel(0);j++) {
        sigma_squared(j) += initial_variance;
        // zero out the initial values
        residuals(j) = 0.0;
    }
    Type persistence = alpha.sum() + beta.sum();
    regressors = v * xi;
    // square y
    vector<Type> residuals_squared = residuals.array().square();
    vector<Type> variance_intercept(timesteps);
    // variance targeting -> not allowed with multiplicative regressors
    Type target_omega = 0.0;
    if (cmodel(3) > 0.5) {
        Type sample_variance = residuals.tail(timesteps - cmodel(0)).square().mean();
        target_omega = sample_variance * (Type(1.0) - persistence);
        // subtract (mean of v) * xi
        vector<Type> meanc = v.colwise().mean();
        Type mean_regressors = (meanc.array() * xi.array()).sum();
        target_omega -= mean_regressors;
        variance_intercept.fill(target_omega);
        ADREPORT(target_omega);
    } else{
        target_omega = omega;
        variance_intercept.fill(omega);
        ADREPORT(target_omega);
    }
    ADREPORT(persistence);
    variance_intercept.array() = variance_intercept.array() + regressors.array();
    // multiplicative regressors
    if (cmodel(4) > 0.5) variance_intercept = variance_intercept.array().exp();
    for(int i = cmodel(0);i<timesteps;i++){
        sigma_squared(i) += variance_intercept(i);
        for(j = 0;j<cmodel(1);j++){
            sigma_squared(i) += alpha(j) * residuals_squared(i - j - 1);
        }
        for(j = 0;j<cmodel(2);j++){
            sigma_squared(i) += beta(j) * sigma_squared(i - j - 1);
        }
    }
    vector<Type> sigma = sigma_squared.sqrt();
    vector<Type> std_residuals = residuals.array() * (Type(1.0)/sigma.array());
    vector<Type> tmp_vector = distfun::distlike(std_residuals, distribution(0), distribution(1), distribution(2), cmodel(5))/sigma.array();
    vector<Type> ll_vector = tmp_vector.tail(timesteps - cmodel(0));
    REPORT(target_omega);
    REPORT(alpha);
    REPORT(beta);
    REPORT(persistence);
    REPORT(initial_variance);
    REPORT(sigma);
    REPORT(ll_vector);
    Type nll = Type(-1.0) * ll_vector.log().sum();
    return(nll);
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
