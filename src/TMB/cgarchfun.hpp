/// @file cgarchfun.hpp
#ifndef cgarchfun_hpp
#define cgarchfun_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type cgarchfun(objective_function<Type>* obj) {
    DATA_VECTOR(y);
    // variance initialization
    DATA_SCALAR(backcast_lambda);
    DATA_INTEGER(samplen);
    DATA_STRING(initmethod);
    PARAMETER(mu);
    PARAMETER(omega);
    PARAMETER(rho);
    PARAMETER(phi);
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
    vector<Type> permanent_component(timesteps);
    vector<Type> transitory_component(timesteps);
    sigma_squared.setZero();
    regressors.setZero();
    permanent_component.setZero();
    transitory_component.setZero();
    int m = v.cols();
    int j = 0;

    // re-scale parameters
    int k = 0;
    mu *= pscale(k);
    k += 1;
    omega *= pscale(k);
    k += 1;
    rho *= pscale(k);
    k += 1;
    phi *= pscale(k);
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

    // variance_targeting (need value of target omega before generating the initialization for the permanent component)
    vector<Type> residuals = y.array() - mu;
    vector<Type> tmp_block = residuals.tail(timesteps - cmodel(0));
    Type initial_variance = garchextra::init_power_variance(tmp_block, initmethod, backcast_lambda, Type(2.0), samplen);
    vector<Type> variance_intercept(timesteps);
    // transitory persistence
    Type persistence = alpha.sum() + beta.sum();
    Type target_omega = 0.0;
    vector<Type> residuals_squared = residuals.array().square();
    // initialize variance
    Type initial_arch = 0.0;
    for(j = 0;j<cmodel(0);j++) {
        // Set to the long run intercept \omega/(1-\rho)
        permanent_component(j) = initial_variance;
        sigma_squared(j) += initial_variance;
        // zero out the initial values
        residuals(j) = 0.0;
        residuals_squared(j) = initial_variance;
        transitory_component(j) = 0.0;
    }

    regressors = v * xi;

    if (cmodel(3) > 0.5) {
        // not using the same variance choice used for initialization as we want a more
        // consistent estimator
        Type sample_variance = residuals.tail(timesteps - cmodel(0)).square().mean();
        target_omega = sample_variance * (Type(1.0) - rho);
        vector<Type> vt_meanc = v.bottomRows(timesteps - cmodel(0)).colwise().mean();
        Type vt_mean_regressors = (vt_meanc.array() * xi.array()).sum();
        // subtract (mean of v) * xi
        target_omega -= vt_mean_regressors;
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
        permanent_component(i) += variance_intercept(i) +  rho * permanent_component(i - 1) + phi * (residuals_squared(i - 1) - sigma_squared(i - 1));
        for(j = 0;j<cmodel(1);j++){
            transitory_component(i) += alpha(j) * (residuals_squared(i - j - 1) - sigma_squared(i - j - 1)) + alpha(j) * transitory_component(i - j - 1);
        }
        for(j = 0;j<cmodel(2);j++){
            transitory_component(i) += beta(j) * transitory_component(i - j - 1);
        }
        sigma_squared(i) += permanent_component(i) + transitory_component(i);
    }
    vector<Type> sigma = sigma_squared.sqrt();
    vector<Type> std_residuals = residuals.array() * (Type(1.0)/sigma.array());
    vector<Type> tmp_vector = distfun::distlike(std_residuals, distribution(0), distribution(1), distribution(2), cmodel(5))/sigma.array();
    vector<Type> ll_vector = tmp_vector.tail(timesteps - cmodel(0));
    REPORT(target_omega);
    REPORT(initial_variance);
    REPORT(initial_arch);
    REPORT(permanent_component);
    REPORT(transitory_component);
    REPORT(sigma);
    REPORT(ll_vector);
    Type nll = Type(-1.0) * ll_vector.log().sum();
    return(nll);
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
