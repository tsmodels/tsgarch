/// @file gjrgarchfun.hpp
#ifndef gjrgarchfun_hpp
#define gjrgarchfun_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type gjrgarchfun(objective_function<Type>* obj) {
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
    vector<Type> negative_indicator(timesteps);
    sigma_squared.setZero();
    sigma.setZero();
    std_residuals.setZero();
    regressors.setZero();
    negative_indicator.setZero();
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

    // variance and arch initialization based on user choice
    vector<Type> residuals = y.array() - mu;
    vector<Type> tmp_block = residuals.tail(timesteps - cmodel(0));
    Type initial_variance = garchextra::init_power_variance(tmp_block, initmethod, backcast_lambda, Type(2.0), samplen);
    Type tmp_arch_init = garchextra::init_gjr(tmp_block, initmethod, backcast_lambda, samplen);
    vector<Type> initial_arch(cmodel(1));
    Type kappa = gjrkappa::gjrgarch_moment_func(distribution(0), distribution(1), distribution(2), dclass);
    ADREPORT(kappa);
    vector<Type> residuals_squared = residuals.array().square();
    for(j = 0;j<cmodel(0);j++) {
        sigma_squared(j) += initial_variance;
        // zero out the initial values
        residuals(j) = 0.0;
        residuals_squared(j) = initial_variance;
        sigma(j) = sqrt(sigma_squared(j));
        negative_indicator(j) = 1.0;
        initial_arch(j) = tmp_arch_init;
    }

    //persistence
    Type persistence = alpha.sum() + beta.sum()  + (gamma.array() * kappa).sum();

    // variance targeting (will not respect user choice of initialization since
    // we use the full sample to capture unconditional sigma)
    regressors = v * xi;
    vector<Type> variance_intercept(timesteps);
    Type target_omega = 0.0;
    if (cmodel(3) > 0.5) {
        Type sample_variance = residuals.tail(timesteps - cmodel(0)).square().mean();
        target_omega = sample_variance * (Type(1.0) - persistence);
        // extract the actual, not zero augmented vector for calculations
        vector<Type> meanc = v.bottomRows(timesteps - cmodel(0)).colwise().mean();
        Type mean_regressors = (meanc.array() * xi.array()).sum();
        // subtract (mean of v) * xi
        target_omega -= mean_regressors;
        variance_intercept.fill(target_omega);
        ADREPORT(target_omega);
    } else {
        target_omega = omega;
        variance_intercept.fill(omega);
        ADREPORT(target_omega);
    }
    ADREPORT(persistence);

    // variance intercept
    variance_intercept.array() = variance_intercept.array() + regressors.array();
    // multiplicative adjustment
    if (cmodel(4) > 0.5) variance_intercept = variance_intercept.array().exp();
    // initialize indicator
    Type indicator = 0.0;

    for(int i = cmodel(0);i<timesteps;i++){
        sigma_squared(i) += variance_intercept(i);
        for(j = 0;j<cmodel(1);j++){
            if((cmodel(1) + j) >= i ) {
                indicator = CppAD::CondExpLe(residuals(i - j - 1), Type(0.0), Type(1.0), Type(0.0));
                negative_indicator(i - j - 1) = indicator;
                sigma_squared(i) += alpha(j) * residuals_squared(i - j - 1) + gamma(j) * initial_arch(j);
            } else {
                indicator = CppAD::CondExpLe(residuals(i - j - 1), Type(0.0), Type(1.0), Type(0.0));
                sigma_squared(i) += alpha(j) * residuals_squared(i - j - 1) + gamma(j) * (indicator * residuals_squared(i - j - 1));
                negative_indicator(i - j - 1) = indicator;
            }
        }
        for(j = 0;j<cmodel(2);j++){
            sigma_squared(i) += beta(j) * sigma_squared(i - j - 1);
        }
        sigma(i) = sqrt(sigma_squared(i));
        std_residuals(i) = residuals(i)/sigma(i);
    }
    vector<Type> tmp_vector = distfun::distlike(std_residuals, distribution(0), distribution(1), distribution(2), dclass)/sigma.array();
    // remove initialization values
    vector<Type> ll_vector = tmp_vector.tail(timesteps - cmodel(0));
    REPORT(target_omega);
    REPORT(negative_indicator);
    REPORT(kappa);
    REPORT(initial_variance);
    REPORT(initial_arch);
    REPORT(sigma);
    REPORT(ll_vector);
    Type nll = Type(-1.0) * ll_vector.log().sum();
    return(nll);
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
