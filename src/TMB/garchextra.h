namespace garchextra {
    template <class Type>
    Type init_power_variance(vector<Type> x, std::string method, Type lambda, Type delta, int samplen)
    {
        const int n = x.rows();
        Type initv = 0.0;
        if (method == "unconditional") {
            initv = x.abs().pow(delta).mean();
        } else if (method == "sample") {
            initv = x.head(samplen).abs().pow(delta).mean();
        } else {
            vector<Type> x_power = x.abs().pow(delta);
            Type mpe = x_power.mean();
            vector<Type> powerlambda(n);
            powerlambda.setZero();
            for(int i = 0;i<=(n-1);i++){
                powerlambda(i) = pow(lambda, i);
            }
            Type tmp = (powerlambda.array() * x_power.array()).sum();
            initv = pow(lambda, n) * mpe + (Type(1.0) - lambda) * tmp;
        }
        return initv;
    }

    template <class Type>
    Type init_log_variance(vector<Type> x, std::string method, Type lambda, int samplen)
    {
        const int n = x.rows();
        Type initv = 0.0;
        if (method == "unconditional") {
            initv = x.square().log().mean();
        } else if (method == "sample") {
            initv = x.head(samplen).square().log().mean();
        } else {
            vector<Type> x_power = x.square().log();
            Type mpe = x_power.mean();
            vector<Type> powerlambda(n);
            powerlambda.setZero();
            for(int i = 0;i<=(n-1);i++){
                powerlambda(i) = pow(lambda, i);
            }
            Type tmp = (powerlambda.array() * x_power.array()).sum();
            initv = pow(lambda, n) * mpe + (Type(1.0) - lambda) * tmp;
        }
        return initv;
    }
}
