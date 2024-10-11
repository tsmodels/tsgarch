namespace egarchkappa {
    // jsu
    template<class Float>
    struct struct_egarch_jsu {
        typedef Float Scalar; // Required by integrate
        Float skew, shape;    // Data
        Float operator() (Float u) {
            Float ans = 0;
            const int give_log = 0;
            ans += fabs(u) * distfun::fwd_jsu_std(u, skew, shape, give_log);
            if (ans == 0) ans = 0;
            using atomic::tiny_ad::isfinite;
            if (!isfinite(ans)) ans = 0;
            return ans;
        }
        Float marginal() {
            using gauss_kronrod::integrate;
            Float ans = integrate(*this, -INFINITY, INFINITY);
            return ans;
        }
    };
    template<class Float>
    Float eval_egarch_jsu(Float skew, Float shape) {
        struct_egarch_jsu<Float> f = {skew, shape};
        return f.marginal();
    }
    TMB_BIND_ATOMIC(egarchjsu, 11, eval_egarch_jsu(x[0], x[1]))
    template<class Type>
    Type jsu_egarch_moment(Type skew, Type shape) {
            vector<Type> args(3); // Last index reserved for derivative order
            args << skew, shape, 0;
            return egarchkappa::egarchjsu(CppAD::vector<Type>(args))[0];
    }

    //sged
    template<class Float>
    struct struct_egarch_sged {
        typedef Float Scalar; // Required by integrate
        Float skew, shape;    // Data
        Float operator() (Float u) {
            Float ans = 0;
            const int give_log = 0;
            ans += fabs(u) * distfun::fwd_sged_std(u, skew, shape, give_log);
            if (ans == 0) ans = 0;
            using atomic::tiny_ad::isfinite;
            if (!isfinite(ans)) ans = 0;
            return ans;
        }
        Float marginal() {
            using gauss_kronrod::integrate;
            Float ans = integrate(*this, -INFINITY, INFINITY);
            return ans;
        }
    };
    template<class Float>
    Float eval_egarch_sged(Float skew, Float shape) {
        struct_egarch_sged<Float> f = {skew, shape};
        return f.marginal();
    }
    TMB_BIND_ATOMIC(egarchsged, 11, eval_egarch_sged(x[0], x[1]))
    template<class Type>
    Type sged_egarch_moment(Type skew, Type shape) {
        vector<Type> args(3); // Last index reserved for derivative order
        args << skew, shape, 0;
        return egarchkappa::egarchsged(CppAD::vector<Type>(args))[0];
    }

    //sstd
    template<class Float>
    struct struct_egarch_sstd {
        typedef Float Scalar; // Required by integrate
        Float skew, shape;    // Data
        Float operator() (Float u) {
            Float ans = 0;
            const int give_log = 0;
            ans += fabs(u) * distfun::fwd_sstd_std(u, skew, shape, give_log);
            if (ans == 0) ans = 0;
            using atomic::tiny_ad::isfinite;
            if (!isfinite(ans)) ans = 0;
            return ans;
        }
        Float marginal() {
            using gauss_kronrod::integrate;
            Float ans = integrate(*this, -INFINITY, INFINITY);
            return ans;
        }
    };
    template<class Float>
    Float eval_egarch_sstd(Float skew, Float shape) {
        struct_egarch_sstd<Float> f = {skew, shape};
        return f.marginal();
    }
    TMB_BIND_ATOMIC(egarchsstd, 11, eval_egarch_sstd(x[0], x[1]))
        template<class Type>
        Type sstd_egarch_moment(Type skew, Type shape) {
            vector<Type> args(3); // Last index reserved for derivative order
            args << skew, shape, 0;
            return egarchkappa::egarchsstd(CppAD::vector<Type>(args))[0];
    }

    //snorm
    template<class Float>
    struct struct_egarch_snorm {
        typedef Float Scalar; // Required by integrate
        Float skew;    // Data
        Float operator() (Float u) {
            Float ans = 0;
            const int give_log = 0;
            ans += fabs(u) * distfun::fwd_snorm_std(u, skew, give_log);
            if (ans == 0) ans = 0;
            using atomic::tiny_ad::isfinite;
            if (!isfinite(ans)) ans = 0;
            return ans;
        }
        Float marginal() {
            using gauss_kronrod::integrate;
            Float ans = integrate(*this, -INFINITY, INFINITY);
            return ans;
        }
    };
    template<class Float>
    Float eval_egarch_snorm(Float skew) {
        struct_egarch_snorm<Float> f = {skew};
        return f.marginal();
    }
    TMB_BIND_ATOMIC(egarchsnorm, 1, eval_egarch_snorm(x[0]))
    template<class Type>
    Type snorm_egarch_moment(Type skew) {
        vector<Type> args(2); // Last index reserved for derivative order
        args << skew, 0;
        return egarchkappa::egarchsnorm(CppAD::vector<Type>(args))[0];
    }

    //nig
    template<class Float>
    struct struct_egarch_nig {
        typedef Float Scalar; // Required by integrate
        Float skew, shape;    // Data
        Float operator() (Float u) {
            Float ans = 0;
            const int give_log = 0;
            ans += fabs(u) * distfun::fwd_nig_std(u, skew, shape, give_log);
            if (ans == 0) ans = 0;
            using atomic::tiny_ad::isfinite;
            if (!isfinite(ans)) ans = 0;
            return ans;
        }
        Float marginal() {
            using gauss_kronrod::integrate;
            Float ans = integrate(*this, -INFINITY, INFINITY);
            return ans;
        }
    };
    template<class Float>
    Float eval_egarch_nig(Float skew, Float shape) {
        struct_egarch_nig<Float> f = {skew, shape};
        return f.marginal();
    }
    TMB_BIND_ATOMIC(egarchnig, 11, eval_egarch_nig(x[0],x[1]))
    template<class Type>
    Type nig_egarch_moment(Type skew, Type shape) {
        vector<Type> args(3); // Last index reserved for derivative order
        args << skew, shape, 0;
        return egarchkappa::egarchnig(CppAD::vector<Type>(args))[0];
    }


    //ghst
    template<class Float>
    struct struct_egarch_ghst {
        typedef Float Scalar; // Required by integrate
        Float skew, shape;    // Data
        Float operator() (Float u) {
            Float ans = 0;
            const int give_log = 0;
            ans += fabs(u) * distfun::fwd_ghst_std(u, skew, shape, give_log);
            if (ans == 0) ans = 0;
            using atomic::tiny_ad::isfinite;
            if (!isfinite(ans)) ans = 0;
            return ans;
        }
        Float marginal() {
            using gauss_kronrod::integrate;
            Float ans = integrate(*this, -INFINITY, INFINITY);
            return ans;
        }
    };
    template<class Float>
    Float eval_egarch_ghst(Float skew, Float shape) {
        struct_egarch_ghst<Float> f = {skew, shape};
        return f.marginal();
    }
    TMB_BIND_ATOMIC(egarchghst, 11, eval_egarch_ghst(x[0],x[1]))
    template<class Type>
    Type ghst_egarch_moment(Type skew, Type shape) {
        vector<Type> args(3); // Last index reserved for derivative order
        args << skew, shape, 0;
        return egarchkappa::egarchghst(CppAD::vector<Type>(args))[0];
    }

    //gh
    template<class Float>
    struct struct_egarch_gh {
        typedef Float Scalar; // Required by integrate
        Float skew, shape, lambda;    // Data
        Float operator() (Float u) {
            Float ans = 0;
            const int give_log = 0;
            ans += fabs(u) * distfun::fwd_gh_std(u, skew, shape, lambda, give_log);
            if (ans == 0) ans = 0;
            using atomic::tiny_ad::isfinite;
            if (!isfinite(ans)) ans = 0;
            return ans;
        }
        Float marginal() {
            using gauss_kronrod::integrate;
            Float ans = integrate(*this, -INFINITY, INFINITY);
            return ans;
        }
    };
    template<class Float>
    Float eval_egarch_gh(Float skew, Float shape, Float lambda) {
        struct_egarch_gh<Float> f = {skew, shape, lambda};
        return f.marginal();
    }
    TMB_BIND_ATOMIC(egarchgh, 111, eval_egarch_gh(x[0],x[1],x[2]))
    template<class Type>
    Type gh_egarch_moment(Type skew, Type shape, Type lambda) {
        vector<Type> args(4); // Last index reserved for derivative order
        args << skew, shape, lambda, 0;
        return egarchkappa::egarchgh(CppAD::vector<Type>(args))[0];
    }

    // closed form solutions for non skewed distributions
    template<class Type>
    Type std_egarch_moment(Type shape){
        Type abs_mom = (Type(2.0)/sqrt(M_PI))*distfun::mygammafn((shape + Type(1.0))/Type(2.0))/distfun::mygammafn(shape/Type(2.0))*sqrt(shape - Type(2.0))/(shape - Type(1.0));
        return abs_mom;
    }

    template<class Type>
    Type ged_egarch_moment(Type shape){
        Type lam = sqrt(Type(1.0)/(pow(Type(2.0),(Type(2.0)/shape))) * distfun::mygammafn(Type(1.0)/shape)/distfun::mygammafn(Type(3.0)/shape));
        Type abs_mom = (pow(Type(2.0),(Type(1.0)/shape)) * lam) * distfun::mygammafn(Type(2.0)/shape)/distfun::mygammafn(Type(1.0)/shape);
        return abs_mom;
    }
    template<class Type>
    Type egarch_moment_func(Type skew, Type shape, Type lambda, int dclass) {
        Type out = 0.0;
        switch(dclass){
        case 1:
            out = sqrt(Type(2.0)/M_PI);
            break;
        case 2:
            out = std_egarch_moment(shape);
            break;
        case 3:
            out = snorm_egarch_moment(skew);
            break;
        case 4:
            out = sstd_egarch_moment(skew, shape);
            break;
        case 5:
            out = ged_egarch_moment(shape);
            break;
        case 6:
            out = sged_egarch_moment(skew, shape);
            break;
        case 7:
            out = nig_egarch_moment(skew, shape);
            break;
        case 8:
            out = gh_egarch_moment(skew, shape, lambda);
            break;
        case 9:
            out = jsu_egarch_moment(skew, shape);
            break;
        case 10:
            out = ghst_egarch_moment(skew, shape);
            break;
        default:
            out = sqrt(Type(2.0)/M_PI);
        }
        return out;
    }
}
