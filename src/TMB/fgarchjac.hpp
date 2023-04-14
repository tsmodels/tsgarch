/// @file fgarchjac.hpp
#ifndef fgarchjac_hpp
#define fgarchjac_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type fgarchjac(objective_function<Type>* obj) {
    PARAMETER_VECTOR(alpha);
    PARAMETER_VECTOR(gamma);
    PARAMETER_VECTOR(eta);
    PARAMETER_VECTOR(beta);
    PARAMETER(delta);
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
    for(j = 0;j<cmodel(1);j++){
        eta(j) *= pscale(j + k);
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
    distribution(0) *= pscale(k);
    distribution(1) *= pscale(k + 1);
    distribution(2) *= pscale(k + 2);
    Type persistence = sum(beta);
    vector<Type> kappa(cmodel(1));
    for(j = 0;j<cmodel(1);j++){
        kappa(j) = fgarchkappa::fgarch_moment_func(gamma(j), eta(j), delta, distribution(0), distribution(1), distribution(2), dclass);
        persistence += alpha(j) * kappa(j);
    }
    return persistence;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
