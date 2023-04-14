/// @file aparchjac.hpp
#ifndef gjrgarchjac_hpp
#define gjrgarchjac_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type gjrgarchjac(objective_function<Type>* obj) {
    PARAMETER_VECTOR(alpha);
    PARAMETER_VECTOR(gamma);
    PARAMETER_VECTOR(beta);
    PARAMETER_VECTOR(distribution);
    // parameter scaling vector
    DATA_IVECTOR(cmodel);
    DATA_VECTOR(pscale);
    int dclass = cmodel(5);
    // re-scale parameters
    int k = 0;
    int j;
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
    distribution(0) *= pscale(k);
    distribution(1) *= pscale(k + 1);
    distribution(2) *= pscale(k + 2);
    Type persistence = sum(beta) + sum(alpha);
    Type kappa = gjrkappa::gjrgarch_moment_func(distribution(0), distribution(1), distribution(2), dclass);
    persistence += (gamma.array() * kappa).sum();
    return persistence;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
