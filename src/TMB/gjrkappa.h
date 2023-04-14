namespace gjrkappa {
    // jsu
    template<class Float>
    struct struct_gjr_jsu {
        typedef Float Scalar; // Required by integrate
        Float skew, shape;    // Data
        Float operator() (Float u) {
            Float ans = 0;
            const int give_log = 0;
            ans += distfun::fwd_jsu_std(u, skew, shape, give_log);
            if (ans == 0) ans = 0;
            using atomic::tiny_ad::isfinite;
            if (!isfinite(ans)) ans = 0;
            return ans;
        }
        Float marginal() {
            using gauss_kronrod::integrate;
            Float ans = integrate(*this, -INFINITY, 0.0);
            return ans;
        }
    };
    template<class Float>
    Float eval_gjr_jsu(Float skew, Float shape) {
        struct_gjr_jsu<Float> f = {skew, shape};
        return f.marginal();
    }
    TMB_BIND_ATOMIC(gjrjsu, 11, eval_gjr_jsu(x[0], x[1]))
    template<class Type>
    Type jsu_gjrgarch_moment(Type skew, Type shape) {
        vector<Type> args(3); // Last index reserved for derivative order
        args << skew, shape, 0;
        return gjrkappa::gjrjsu(CppAD::vector<Type>(args))[0];
    }

    //sged
    template<class Float>
    struct struct_gjr_sged {
        typedef Float Scalar; // Required by integrate
        Float skew, shape;    // Data
        Float operator() (Float u) {
            Float ans = 0;
            const int give_log = 0;
            ans += distfun::fwd_sged_std(u, skew, shape, give_log);
            if (ans == 0) ans = 0;
            using atomic::tiny_ad::isfinite;
            if (!isfinite(ans)) ans = 0;
            return ans;
        }
        Float marginal() {
            using gauss_kronrod::integrate;
            Float ans = integrate(*this, -INFINITY, 0.0);
            return ans;
        }
    };
    template<class Float>
    Float eval_gjr_sged(Float skew, Float shape) {
        struct_gjr_sged<Float> f = {skew, shape};
        return f.marginal();
    }
    TMB_BIND_ATOMIC(gjrsged, 11, eval_gjr_sged(x[0], x[1]))
    template<class Type>
    Type sged_gjrgarch_moment(Type skew, Type shape) {
        vector<Type> args(3); // Last index reserved for derivative order
        args << skew, shape, 0;
        return gjrkappa::gjrsged(CppAD::vector<Type>(args))[0];
    }

    //sstd
    template<class Float>
    struct struct_gjr_sstd {
        typedef Float Scalar; // Required by integrate
        Float skew, shape;    // Data
        Float operator() (Float u) {
            Float ans = 0;
            const int give_log = 0;
            ans += distfun::fwd_sstd_std(u, skew, shape, give_log);
            if (ans == 0) ans = 0;
            using atomic::tiny_ad::isfinite;
            if (!isfinite(ans)) ans = 0;
            return ans;
        }
        Float marginal() {
            using gauss_kronrod::integrate;
            Float ans = integrate(*this, -INFINITY, 0.0);
            return ans;
        }
    };
    template<class Float>
    Float eval_gjr_sstd(Float skew, Float shape) {
        struct_gjr_sstd<Float> f = {skew, shape};
        return f.marginal();
    }
    TMB_BIND_ATOMIC(gjrsstd, 11, eval_gjr_sstd(x[0], x[1]))
    template<class Type>
    Type sstd_gjrgarch_moment(Type skew, Type shape) {
        vector<Type> args(3); // Last index reserved for derivative order
        args << skew, shape, 0;
        return gjrkappa::gjrsstd(CppAD::vector<Type>(args))[0];
    }

    //snorm
    template<class Float>
    struct struct_gjr_snorm {
        typedef Float Scalar; // Required by integrate
        Float skew;    // Data
        Float operator() (Float u) {
            Float ans = 0;
            const int give_log = 0;
            ans += distfun::fwd_snorm_std(u, skew, give_log);
            if (ans == 0) ans = 0;
            using atomic::tiny_ad::isfinite;
            if (!isfinite(ans)) ans = 0;
            return ans;
        }
        Float marginal() {
            using gauss_kronrod::integrate;
            Float ans = integrate(*this, -INFINITY, 0.0);
            return ans;
        }
    };
    template<class Float>
    Float eval_gjr_snorm(Float skew) {
        struct_gjr_snorm<Float> f = {skew};
        return f.marginal();
    }
    TMB_BIND_ATOMIC(gjrsnorm, 1, eval_gjr_snorm(x[0]))
    template<class Type>
    Type snorm_gjrgarch_moment(Type skew) {
        vector<Type> args(2); // Last index reserved for derivative order
        args << skew, 0;
        return gjrkappa::gjrsnorm(CppAD::vector<Type>(args))[0];
    }

    //nig
    template<class Float>
    struct struct_gjr_nig {
        typedef Float Scalar; // Required by integrate
        Float skew, shape;    // Data
        Float operator() (Float u) {
            Float ans = 0;
            const int give_log = 0;
            ans += distfun::fwd_nig_std(u, skew, shape, give_log);
            if (ans == 0) ans = 0;
            using atomic::tiny_ad::isfinite;
            if (!isfinite(ans)) ans = 0;
            return ans;
        }
        Float marginal() {
            using gauss_kronrod::integrate;
            Float ans = integrate(*this, -INFINITY, 0.0);
            return ans;
        }
    };
    template<class Float>
    Float eval_gjr_nig(Float skew, Float shape) {
        struct_gjr_nig<Float> f = {skew, shape};
        return f.marginal();
    }
    TMB_BIND_ATOMIC(gjrnig, 11, eval_gjr_nig(x[0],x[1]))
    template<class Type>
    Type nig_gjrgarch_moment(Type skew, Type shape) {
        vector<Type> args(3); // Last index reserved for derivative order
        args << skew, shape, 0;
        return gjrkappa::gjrnig(CppAD::vector<Type>(args))[0];
    }


    //ghst
    template<class Float>
    struct struct_gjr_ghst {
        typedef Float Scalar; // Required by integrate
        Float skew, shape;    // Data
        Float operator() (Float u) {
            Float ans = 0;
            const int give_log = 0;
            ans += distfun::fwd_ghst_std(u, skew, shape, give_log);
            if (ans == 0) ans = 0;
            using atomic::tiny_ad::isfinite;
            if (!isfinite(ans)) ans = 0;
            return ans;
        }
        Float marginal() {
            using gauss_kronrod::integrate;
            Float ans = integrate(*this, -INFINITY, 0.0);
            return ans;
        }
    };
    template<class Float>
    Float eval_gjr_ghst(Float skew, Float shape) {
        struct_gjr_ghst<Float> f = {skew, shape};
        return f.marginal();
    }
    TMB_BIND_ATOMIC(gjrghst, 11, eval_gjr_ghst(x[0],x[1]))
    template<class Type>
    Type ghst_gjrgarch_moment(Type skew, Type shape) {
        vector<Type> args(3); // Last index reserved for derivative order
        args << skew, shape, 0;
        return gjrkappa::gjrghst(CppAD::vector<Type>(args))[0];
    }

    //gh
    template<class Float>
    struct struct_gjr_gh {
        typedef Float Scalar; // Required by integrate
        Float skew, shape, lambda;    // Data
        Float operator() (Float u) {
            Float ans = 0;
            const int give_log = 0;
            ans += distfun::fwd_gh_std(u, skew, shape, lambda, give_log);
            if (ans == 0) ans = 0;
            using atomic::tiny_ad::isfinite;
            if (!isfinite(ans)) ans = 0;
            return ans;
        }
        Float marginal() {
            using gauss_kronrod::integrate;
            Float ans = integrate(*this, -INFINITY, 0.0);
            return ans;
        }
    };
    template<class Float>
    Float eval_gjr_gh(Float skew, Float shape, Float lambda) {
        struct_gjr_gh<Float> f = {skew, shape, lambda};
        return f.marginal();
    }
    TMB_BIND_ATOMIC(gjrgh, 111, eval_gjr_gh(x[0],x[1],x[2]))
    template<class Type>
    Type gh_gjrgarch_moment(Type skew, Type shape, Type lambda) {
        vector<Type> args(4); // Last index reserved for derivative order
        args << skew, shape, lambda, 0;
        return gjrkappa::gjrgh(CppAD::vector<Type>(args))[0];
    }

    template<class Type>
    Type gjrgarch_moment_func(Type skew, Type shape, Type lambda, int dclass) {
        Type out = 0.0;
        switch(dclass){
        case 1:
            out = Type(0.5);
            break;
        case 2:
            out = Type(0.5);
            break;
        case 3:
            out = snorm_gjrgarch_moment(skew);
            break;
        case 4:
            out = sstd_gjrgarch_moment(skew, shape);
            break;
        case 5:
            out = Type(0.5);
            break;
        case 6:
            out = sged_gjrgarch_moment(skew, shape);
            break;
        case 7:
            out = nig_gjrgarch_moment(skew, shape);
            break;
        case 8:
            out = gh_gjrgarch_moment(skew, shape, lambda);
            break;
        case 9:
            out = jsu_gjrgarch_moment(skew, shape);
            break;
        case 10:
            out = ghst_gjrgarch_moment(skew, shape);
            break;
        default:
            out = Type(0.5);
        }
        return out;
    }
}
