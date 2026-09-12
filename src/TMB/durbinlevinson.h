/// @file durbinlevinson.h
/// Durbin-Levinson (Jones, 1980) forward transform used to reparameterize the
/// ARMA mean equation. A vector of raw parameters r_1,...,r_n each confined to
/// the open interval (-1,1) maps, via the recursion below, onto AR
/// coefficients guaranteed to lie in the stationarity region (equivalently,
/// applied to the negated raw parameters and negated back, onto MA
/// coefficients guaranteed to lie in the invertibility region). Because the
/// recursion is built entirely from ordinary arithmetic (no branching, no
/// eigendecomposition), it differentiates exactly under CppAD/TMB's automatic
/// differentiation with no need for RTMB or nonlinear inequality constraints.
/// See R/arma.R for the (identical) R-level implementation used to generate
/// starting values.
#ifndef durbinlevinson_h
#define durbinlevinson_h

namespace garchextra {

template<class Type>
vector<Type> durbin_levinson(const vector<Type>& r) {
    int n = r.size();
    vector<Type> phi(n);
    if (n == 0) return phi;
    phi(0) = r(0);
    for (int k = 1; k < n; k++) {
        vector<Type> phi_new(k + 1);
        phi_new(k) = r(k);
        for (int i = 0; i < k; i++) {
            phi_new(i) = phi(i) - r(k) * phi(k - 1 - i);
        }
        phi = phi_new;
    }
    return phi;
}

// AR coefficients from raw (pacf-space) parameters: guaranteed stationary.
template<class Type>
vector<Type> pacf_to_ar(const vector<Type>& r) {
    return durbin_levinson(r);
}

// MA coefficients (y_t = ... + eps_t + theta_1*eps_{t-1} + ...) from raw
// (pacf-space) parameters: guaranteed invertible, via the AR/MA duality of
// the Durbin-Levinson recursion (apply to the negated input, negate output).
template<class Type>
vector<Type> pacf_to_ma(const vector<Type>& r) {
    vector<Type> neg_r = -r;
    vector<Type> theta = durbin_levinson(neg_r);
    return (Type(-1.0) * theta);
}

} // namespace garchextra

#endif
