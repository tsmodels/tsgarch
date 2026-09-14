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
    // raw (pacf-space) Durbin-Levinson parameters for the ARMA mean equation;
    // see durbinlevinson.h. Guaranteed stationary/invertible for any value in
    // (-1,1), so no additional nonlinear constraint is required.
    PARAMETER_VECTOR(arpacf);
    PARAMETER_VECTOR(mapacf);
    // mean-equation regressor coefficients
    PARAMETER_VECTOR(tau);
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
    // mean equation regressors
    DATA_MATRIX(x);
    // model flags [maxpq arch_order garch_order variance_targeting multiplicative distribution_no ar_order ma_order armax]
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
    // tau rows always number exactly x.cols() (>=1: the dummy column/row
    // convention matches the arpacf/mapacf dummy convention), so no
    // zero-column special case is needed here
    const int mx = x.cols();
    for(j = 0;j<mx;j++) { tau(j) *= pscale(j + k); }
    k += mx;
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

    // ARMA mean equation:
    //   (y_t - mu) = sum_i phi_i * (y_{t-i} - mu) + eps_t + sum_j theta_j * eps_{t-j}
    // so that mu retains its interpretation as the unconditional mean of y
    // (when ar_order = ma_order = 0 this reduces exactly to eps_t = y_t - mu,
    // i.e. the pre-ARMA behavior). arma_ar/arma_ma are guaranteed stationary/invertible.
    vector<Type> arma_ar = garchextra::pacf_to_ar(arpacf);
    vector<Type> arma_ma = garchextra::pacf_to_ma(mapacf);
    vector<Type> xtau = x * tau;
    const int armax = cmodel(8);
    vector<Type> z(timesteps);
    if (armax > 0) {
        z = y.array() - mu;
    } else {
        z = y.array() - mu - xtau.array();
    }
    // y's and x's pre-sample rows (indices < cmodel(0)) are zero-padded by
    // the R wrapper, which would otherwise leak z(presample) = 0 - mu = -mu
    // into the AR feedback for the first ar_order real observations (and
    // xtau = 0 there in any case). Force the pre-sample z (and residuals,
    // already zero below) to 0, i.e. "the process starts at its
    // unconditional mean with zero shocks".
    for (int i = 0; i < cmodel(0); i++) z(i) = Type(0.0);
    vector<Type> residuals(timesteps);
    residuals.setZero();
    // conditional_mean(i) is the model's fitted conditional mean of y at
    // time i (i.e. mu + the regressor contribution + the AR/MA deviation
    // term), so that
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
        if (armax > 0) mean_i += xtau(i);
        residuals(i) = z(i) - mean_i;
        // y(i) - residuals(i) is algebraically identical to mu + mean_i
        // under arma_errors and is the correct conditional mean under
        // armax (where xtau enters the mean recursion itself)
        conditional_mean(i) = y(i) - residuals(i);
    }
    // variance and arch initialization based on user choice
    vector<Type> tmp_block = residuals.tail(timesteps - cmodel(0));
    Type initial_variance = garchextra::init_power_variance(tmp_block, initmethod, backcast_lambda, Type(2.0), samplen);
    Type initial_log_variance = log(initial_variance);
    vector<Type> initial_arch(cmodel(1));
    for(j = 0;j<cmodel(0);j++) {
        sigma_squared(j) += initial_variance;
        log_sigma_squared(j) = initial_log_variance;
        // zero out the initial values
        residuals(j) = 0.0;
        sigma(j) = sqrt(sigma_squared(j));
    }
    // initial_arch is only ever indexed up to cmodel(1) (the ARCH order),
    // which can now be smaller than cmodel(0) (the combined GARCH/ARMA
    // pre-sample length) - initialize it in its own, correctly-bounded loop
    // rather than the cmodel(0)-bounded loop above.
    for(j = 0;j<cmodel(1);j++) {
        initial_arch(j) = 0.0;
    }

    // expectation of abs(z)
    Type kappa = egarchkappa::egarch_moment_func(distribution(0), distribution(1), distribution(2), dclass);
    ADREPORT(kappa);

    // persistence
    Type persistence = beta.sum();

    // variance targeting (will not respect user choice of initialization since
    // we use the full sample to capture unconditional sigma)
    regressors = v * xi;
    vector<Type> residuals_squared = residuals.array().square();
    vector<Type> variance_intercept(timesteps);
    Type target_omega = 0.0;
    if (cmodel(3) > 0.5) {
        Type log_sample_variance = log(residuals.tail(timesteps - cmodel(0)).square().mean());
        target_omega = log_sample_variance * (Type(1.0) - persistence);
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

    for(int i = cmodel(0);i<timesteps;i++){
        log_sigma_squared(i) += variance_intercept(i);
        for(j = 0;j<cmodel(1);j++){
            if ((cmodel(1) + j) >= i) {
                log_sigma_squared(i) += alpha(j) * Type(0.0) + gamma(j) * initial_arch(j);
            } else {
                log_sigma_squared(i) += alpha(j) * std_residuals(i - j - 1) + gamma(j) * (fabs(std_residuals(i - j - 1)) - kappa);
            }
        }
        for(j = 0;j<cmodel(2);j++){
            log_sigma_squared(i) += beta(j) * log_sigma_squared(i - j - 1);
        }
        sigma_squared(i) = exp(log_sigma_squared(i));
        sigma(i) = sqrt(sigma_squared(i));
        std_residuals(i) = residuals(i)/sigma(i);
    }
    vector<Type> tmp_vector = distfun::distlike(std_residuals, distribution(0), distribution(1), distribution(2), dclass)/sigma.array();
    // remove initialization values
    vector<Type> ll_vector = tmp_vector.tail(timesteps - cmodel(0));
    REPORT(target_omega);
    REPORT(kappa);
    REPORT(initial_variance);
    REPORT(initial_arch);
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
