namespace garchextra {
    // Eigen 5.0 removed the free pow(ArrayBase, scalar) for non-arithmetic
    // (TMB CppAD::AD<double>) scalar exponents. Build the same elementwise
    // pow-by-scalar via unaryExpr, which works across Eigen 3.4 and 5.0.
    template <class Derived, class ExponentScalar>
    inline auto pow_cwise(const Eigen::ArrayBase<Derived>& x, const ExponentScalar& e) {
        typedef typename Derived::Scalar Scalar;
        return x.derived().unaryExpr([e](const Scalar& v) { return Eigen::numext::pow(v, e); });
    }

    template <class Type>
    Type init_power_variance(vector<Type> x, std::string method, Type lambda, Type delta, int samplen)
    {
        const int n = x.rows();
        Type initv = 0.0;
        if (method == "unconditional") {
            initv = x.abs().pow(2.0).mean();
        } else if (method == "sample") {
            initv = x.head(samplen).abs().pow(2.0).mean();
        } else {
            vector<Type> x_power = x.abs().pow(2.0);
            Type mpe = x_power.mean();
            vector<Type> powerlambda(n);
            powerlambda.setZero();
            for(int i = 0;i<=(n-1);i++){
                powerlambda(i) = pow(lambda, i);
            }
            Type tmp = (powerlambda.array() * x_power.array()).sum();
            initv = pow(lambda, n) * mpe + (Type(1.0) - lambda) * tmp;
        }
        initv = pow(initv, delta/2.0);
        return initv;
    }

    template <class Type>
    Type init_aparch(vector<Type> x, std::string method, Type gamma, Type delta, Type lambda, int samplen)
    {
        const int n = x.rows();
        Type inite = 0.0;
        if (method == "unconditional") {
            inite += pow_cwise(x.abs().array() - gamma * x.array(), delta).mean();
        } else if (method == "sample") {
            inite += pow_cwise(x.head(samplen).abs().array() - gamma * x.head(samplen).array(), delta).mean();
        } else {
            vector<Type> x_power = pow_cwise(x.abs().array() - gamma * x.array(), delta);
            Type mean_x_power = x_power.mean();
            vector<Type> powerlambda(n);
            powerlambda.setZero();
            for(int i = 0;i<=(n-1);i++){
                powerlambda(i) = pow(lambda, i);
            }
            Type tmp = (powerlambda.array() * x_power.array()).sum();
            inite = pow(lambda, n) * mean_x_power + (Type(1.0) - lambda) * tmp;
        }
        return inite;
    }

    template <class Type>
    Type init_fgarch(vector<Type> x, std::string method, Type gamma, Type eta, Type delta, Type lambda, int samplen)
    {
        const int n = x.rows();
        Type inite = 0.0;
        if (method == "unconditional") {
            inite += pow_cwise((x.array() - eta).abs() - gamma * (x.array() - eta), delta).mean();
        } else if (method == "sample") {
            inite += pow_cwise((x.head(samplen).array() - eta).abs() - gamma * (x.head(samplen).array() - eta), delta).mean();
        } else {
            vector<Type> x_power = pow_cwise((x.array() - eta).abs() - gamma * (x.array() - eta), delta);
            Type mean_x_power = x_power.mean();
            vector<Type> powerlambda(n);
            powerlambda.setZero();
            for(int i = 0;i<=(n-1);i++){
                powerlambda(i) = pow(lambda, i);
            }
            Type tmp = (powerlambda.array() * x_power.array()).sum();
            inite = pow(lambda, n) * mean_x_power + (Type(1.0) - lambda) * tmp;
        }
        return inite;
    }

    template <class Type>
    Type init_gjr(vector<Type> x, std::string method, Type lambda, int samplen)
    {
        const int n = x.rows();
        Type inite = 0.0;
        if (method == "unconditional") {
            vector<Type> vsign = (Type(1.0) - x.array().sign())/Type(2.0);
            inite += (x.pow(2.0).array() * vsign.array()).mean();
        } else if (method == "sample") {
            vector<Type> vsign = (Type(1.0) - x.head(samplen).array().sign())/Type(2.0);
            inite += (x.head(samplen).pow(2.0).array() * vsign.array()).mean();
        } else {
            vector<Type> vsign = (Type(1.0) - x.array().sign())/Type(2.0);
            vector<Type> nege  = (x.pow(2.0).array() * vsign.array());
            Type mean_nege = nege.mean();
            vector<Type> powerlambda(n);
            powerlambda.setZero();
            for(int i = 0;i<=(n-1);i++){
                powerlambda(i) = pow(lambda, i);
            }
            Type tmp = (powerlambda.array() * nege.array()).sum();
            inite = pow(lambda, n) * mean_nege + (Type(1.0) - lambda) * tmp;
        }
        return inite;
    }

    template <class Type>
    vector<Type> init_mean(matrix<Type> x, std::string method, int samplen)
    {
        const int m = x.cols();
        vector<Type> x_mean(m);
        if (method == "unconditional") {
            x_mean = x.colwise().mean();
        } else if (method == "sample") {
            x_mean = x.topRows(samplen).colwise().mean();
        } else {
            x_mean = x.colwise().mean();
        }
        return x_mean;
    }
}
