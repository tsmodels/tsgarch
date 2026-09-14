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
    PARAMETER(delta);
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
    delta *= pscale(k);
    k += 1;
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
    // variance initialization based on user choice
    // extract the actual, not zero augmented vector for calculations
    vector<Type> tmp_block = residuals.tail(timesteps - cmodel(0));
    Type initial_power_sigma = garchextra::init_power_variance(tmp_block, initmethod, backcast_lambda, delta, samplen);
    Type initial_variance = pow(initial_power_sigma, Type(2.0)/delta);
    for(j = 0;j<cmodel(0);j++) {
        sigma_power(j) += initial_power_sigma;
        // zero out the initial values
        residuals(j) = 0.0;
        sigma(j) = pow(sigma_power(j), Type(1.0)/delta);
    }

    // arch initialization (respects user choice for initialization)
    // persistence
    Type persistence = beta.sum();
    vector<Type> kappa(cmodel(1));
    vector<Type> initial_arch(cmodel(1));
    initial_arch.setZero();
    for(j = 0;j<cmodel(1);j++) {
        kappa(j) = aparchkappa::aparch_moment_func(gamma(j), delta, distribution(0), distribution(1), distribution(2), dclass);
        persistence += alpha(j) * kappa(j);
        initial_arch(j) = garchextra::init_aparch(tmp_block, initmethod, gamma(j), delta, backcast_lambda, samplen);
    }
    ADREPORT(persistence);
    ADREPORT(kappa);

    // variance targeting (will not respect user choice of initialization since
    // we use the full sample to capture unconditional sigma)
    regressors = v * xi;
    vector<Type> power_sigma_intercept(timesteps);
    Type target_omega = 0.0;
    if (cmodel(3) > 0.5) {
        Type sample_power_sigma = residuals.tail(timesteps - cmodel(0)).pow(2.0).mean();
        sample_power_sigma = pow(sample_power_sigma, delta/2.0);
        target_omega = sample_power_sigma * (Type(1.0) - persistence);
        // use the actual, not zero augmented vector for calculations
        vector<Type> meanc = v.bottomRows(timesteps - cmodel(0)).colwise().mean();
        Type mean_regressors = (meanc.array() * xi.array()).sum();
        // subtract (mean of v) * xi
        target_omega -= mean_regressors;
        power_sigma_intercept.fill(target_omega);
        ADREPORT(target_omega);
    } else {
        target_omega = omega;
        power_sigma_intercept.fill(omega);
        ADREPORT(target_omega);
    }

    // variance intercept
    power_sigma_intercept.array() = power_sigma_intercept.array() + regressors.array();
    // multiplicative adjustment
    if (cmodel(4) > 0.5) power_sigma_intercept = power_sigma_intercept.array().exp();
    for(int i = cmodel(0);i<timesteps;i++){
        sigma_power(i) += power_sigma_intercept(i);
        for(j = 0;j<cmodel(1);j++){
            if((cmodel(1) + j) >= i ) {
                sigma_power(i) += alpha(j) * initial_arch(j);
            } else {
                sigma_power(i) += alpha(j) * pow(fabs(residuals(i - j - 1)) - gamma(j) * residuals(i - j - 1), delta);
            }
        }
        for(j = 0;j<cmodel(2);j++){
            sigma_power(i) += beta(j) * sigma_power(i - j - 1);
        }
        sigma(i) = pow(sigma_power(i), Type(1.0)/delta);
        std_residuals(i) = residuals(i)/sigma(i);
    }
    vector<Type> tmp_vector = distfun::distlike(std_residuals, distribution(0), distribution(1), distribution(2), dclass)/sigma.array();
    // remove initialization values
    vector<Type> ll_vector = tmp_vector.tail(timesteps - cmodel(0));
    REPORT(target_omega);
    REPORT(persistence);
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
