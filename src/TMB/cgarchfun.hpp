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
    // raw (pacf-space) Durbin-Levinson parameters for the ARMA mean equation;
    // see durbinlevinson.h. Guaranteed stationary/invertible for any value in
    // (-1,1), so no additional nonlinear constraint is required.
    PARAMETER_VECTOR(arpacf);
    PARAMETER_VECTOR(mapacf);
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
    // model flags [maxpq arch_order garch_order variance_targeting multiplicative distribution_no ar_order ma_order]
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
    const int ar_order = cmodel(6);
    const int ma_order = cmodel(7);

    // re-scale parameters
    int k = 0;
    mu *= pscale(k);
    k += 1;
    for(j = 0;j<ar_order;j++) {
        arpacf(j) *= pscale(j + k);
    }
    if (ar_order == 0) {
        k += 1;
    } else {
        k += ar_order;
    }
    for(j = 0;j<ma_order;j++) {
        mapacf(j) *= pscale(j + k);
    }
    if (ma_order == 0) {
        k += 1;
    } else {
        k += ma_order;
    }
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

    // ARMA mean equation:
    //   (y_t - mu) = sum_i phi_i * (y_{t-i} - mu) + eps_t + sum_j theta_j * eps_{t-j}
    // so that mu retains its interpretation as the unconditional mean of y
    // (when ar_order = ma_order = 0 this reduces exactly to eps_t = y_t - mu,
    // i.e. the pre-ARMA behavior). arma_ar/arma_ma are guaranteed stationary/invertible.
    // (named arma_ar/arma_ma rather than phi/theta since this model already
    // has its own scalar PARAMETER(phi) for the transitory-component coefficient)
    vector<Type> arma_ar = garchextra::pacf_to_ar(arpacf);
    vector<Type> arma_ma = garchextra::pacf_to_ma(mapacf);
    vector<Type> z = y.array() - mu;
    // y's pre-sample rows (indices < cmodel(0)) are zero-padded by the R
    // wrapper, which would otherwise leak z(presample) = 0 - mu = -mu into
    // the AR feedback for the first ar_order real observations. Force the
    // pre-sample z (and residuals, already zero below) to 0, i.e. "the
    // process starts at its unconditional mean with zero shocks".
    for (int i = 0; i < cmodel(0); i++) z(i) = Type(0.0);
    vector<Type> residuals(timesteps);
    residuals.setZero();
    // conditional_mean(i) is the model's fitted conditional mean of y at
    // time i (i.e. mu + the AR/MA deviation term), so that
    // residuals(i) = y(i) - conditional_mean(i) always holds exactly - the
    // same relationship used by fitted()/residuals() on the R side. Reduces
    // to a constant mu everywhere when ar_order = ma_order = 0. Pre-sample
    // rows are set to mu (consistent with the "start at the unconditional
    // mean" convention used for the pre-sample z/residuals above).
    vector<Type> conditional_mean(timesteps);
    conditional_mean.fill(mu);
    for (int i = cmodel(0); i < timesteps; i++) {
        Type mean_i = Type(0.0);
        for (j = 0; j < ar_order; j++) {
            mean_i += arma_ar(j) * z(i - j - 1);
        }
        for (j = 0; j < ma_order; j++) {
            mean_i += arma_ma(j) * residuals(i - j - 1);
        }
        residuals(i) = z(i) - mean_i;
        conditional_mean(i) = mu + mean_i;
    }
    // variance_targeting (need value of target omega before generating the initialization for the permanent component)
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
    REPORT(arma_ar);
    REPORT(arma_ma);
    REPORT(residuals);
    REPORT(conditional_mean);
    ADREPORT(arma_ar);
    ADREPORT(arma_ma);
    Type nll = Type(-1.0) * ll_vector.log().sum();
    return(nll);
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
