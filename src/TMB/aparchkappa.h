namespace aparchkappa {
    // jsu
    template<class Float>
    struct struct_aparch_jsu {
        typedef Float Scalar; // Required by integrate
        Float gamma, delta, skew, shape;    // Data
        Float operator() (Float u) {
            Float ans = 0;
            const int give_log = 0;
            ans += pow(fabs(u) - gamma * u, delta) * distfun::fwd_jsu_std(u, skew, shape, give_log);
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
    Float eval_aparch_jsu(Float gamma, Float delta, Float skew, Float shape) {
        struct_aparch_jsu<Float> f = {gamma, delta, skew, shape};
        return f.marginal();
    }

    TMB_BIND_ATOMIC(aparchjsu, 1111, eval_aparch_jsu(x[0], x[1], x[2], x[3]))

    template<class Type>
    Type jsu_aparch_moment(Type gamma, Type delta, Type skew, Type shape) {
            vector<Type> args(5); // Last index reserved for derivative order
            args << gamma, delta, skew, shape, 0;
            return aparchkappa::aparchjsu(CppAD::vector<Type>(args))[0];
    }

    //sged
    template<class Float>
    struct struct_aparch_sged {
        typedef Float Scalar; // Required by integrate
        Float gamma, delta, skew, shape;    // Data
        Float operator() (Float u) {
            Float ans = 0;
            const int give_log = 0;
            ans += pow(fabs(u) - gamma * u, delta) * distfun::fwd_sged_std(u, skew, shape, give_log);
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
    Float eval_aparch_sged(Float gamma, Float delta, Float skew, Float shape) {
        struct_aparch_sged<Float> f = {gamma, delta, skew, shape};
        return f.marginal();
    }

    TMB_BIND_ATOMIC(aparchsged, 1111, eval_aparch_sged(x[0], x[1], x[2], x[3]))

    template<class Type>
    Type sged_aparch_moment(Type gamma, Type delta, Type skew, Type shape) {
        vector<Type> args(5); // Last index reserved for derivative order
        args << gamma, delta, skew, shape, 0;
        return aparchkappa::aparchsged(CppAD::vector<Type>(args))[0];
    }

    //sstd
    template<class Float>
    struct struct_aparch_sstd {
        typedef Float Scalar; // Required by integrate
        Float gamma, delta, skew, shape;    // Data
        Float operator() (Float u) {
            Float ans = 0;
            const int give_log = 0;
            ans += pow(fabs(u) - gamma * u, delta) * distfun::fwd_sstd_std(u, skew, shape, give_log);
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
    Float eval_aparch_sstd(Float gamma, Float delta, Float skew, Float shape) {
        struct_aparch_sstd<Float> f = {gamma, delta, skew, shape};
        return f.marginal();
    }

    TMB_BIND_ATOMIC(aparchsstd, 1111, eval_aparch_sstd(x[0], x[1], x[2], x[3]))

    template<class Type>
    Type sstd_aparch_moment(Type gamma, Type delta, Type skew, Type shape) {
        vector<Type> args(5); // Last index reserved for derivative order
        args << gamma, delta, skew, shape, 0;
        return aparchkappa::aparchsstd(CppAD::vector<Type>(args))[0];
    }

    //snorm
    template<class Float>
    struct struct_aparch_snorm {
        typedef Float Scalar; // Required by integrate
        Float gamma, delta, skew;    // Data
        Float operator() (Float u) {
            Float ans = 0;
            const int give_log = 0;
            ans += pow(fabs(u) - gamma * u, delta) * distfun::fwd_snorm_std(u, skew, give_log);
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
    Float eval_aparch_snorm(Float gamma, Float delta, Float skew) {
        struct_aparch_snorm<Float> f = {gamma, delta, skew};
        return f.marginal();
    }
    TMB_BIND_ATOMIC(aparchsnorm, 111, eval_aparch_snorm(x[0], x[1], x[2]))
    template<class Type>
    Type snorm_aparch_moment(Type gamma, Type delta, Type skew) {
        vector<Type> args(4); // Last index reserved for derivative order
        args << gamma, delta, skew, 0;
        return aparchkappa::aparchsnorm(CppAD::vector<Type>(args))[0];
    }

    //nig
    template<class Float>
    struct struct_aparch_nig {
        typedef Float Scalar; // Required by integrate
        Float gamma, delta, skew, shape;    // Data
        Float operator() (Float u) {
            Float ans = 0;
            const int give_log = 0;
            ans += pow(fabs(u) - gamma * u, delta) * distfun::fwd_nig_std(u, skew, shape, give_log);
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
    Float eval_aparch_nig(Float gamma, Float delta, Float skew, Float shape) {
        struct_aparch_nig<Float> f = {gamma, delta, skew, shape};
        return f.marginal();
    }

    TMB_BIND_ATOMIC(aparchnig, 1111, eval_aparch_nig(x[0], x[1], x[2], x[3]))

    template<class Type>
    Type nig_aparch_moment(Type gamma, Type delta, Type skew, Type shape) {
        vector<Type> args(5); // Last index reserved for derivative order
        args << gamma, delta, skew, shape, 0;
        return aparchkappa::aparchnig(CppAD::vector<Type>(args))[0];
    }


      template<class Float>
      struct struct_aparch_ghst {
          typedef Float Scalar; // Required by integrate
          Float gamma, delta, skew, shape;    // Data
          Float operator() (Float u) {
              Float ans = 0;
              const int give_log = 0;
              ans += pow(fabs(u) - gamma * u, delta) * distfun::fwd_ghst_std(u, skew, shape, give_log);
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
      Float eval_aparch_ghst(Float gamma, Float delta, Float skew, Float shape) {
          struct_aparch_ghst<Float> f = {gamma, delta, skew, shape};
          return f.marginal();
      }

      TMB_BIND_ATOMIC(aparchghst, 1111, eval_aparch_ghst(x[0], x[1], x[2], x[3]))

      template<class Type>
      Type ghst_aparch_moment(Type gamma, Type delta, Type skew, Type shape) {
          vector<Type> args(5); // Last index reserved for derivative order
          args << gamma, delta, skew, shape, 0;
          return aparchkappa::aparchghst(CppAD::vector<Type>(args))[0];
      }

      template<class Float>
      struct struct_aparch_gh {
          typedef Float Scalar; // Required by integrate
          Float gamma, delta, skew, shape, lambda;    // Data
          Float operator() (Float u) {
              Float ans = 0;
              const int give_log = 0;
              ans += pow(fabs(u) - gamma * u, delta) * distfun::fwd_gh_std(u, skew, shape, lambda, give_log);
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
      Float eval_aparch_gh(Float gamma, Float delta, Float skew, Float shape, Float lambda) {
          struct_aparch_gh<Float> f = {gamma, delta, skew, shape, lambda};
          return f.marginal();
      }

      TMB_BIND_ATOMIC(aparchgh, 11111, eval_aparch_gh(x[0], x[1], x[2], x[3], x[4]))

      template<class Type>
      Type gh_aparch_moment(Type gamma, Type delta, Type skew, Type shape, Type lambda) {
          vector<Type> args(6); // Last index reserved for derivative order
          args << gamma, delta, skew, shape, lambda, 0;
          return aparchkappa::aparchgh(CppAD::vector<Type>(args))[0];
      }

    // we have closed form solutions for the non skewed distributions
    template<class Type>
    Type norm_aparch_moment(Type gamma, Type delta){
        Type bc_mom = (Type(1.0)/sqrt(M_PI)) * (pow(Type(1.0) - gamma, delta) + pow(Type(1.0) + gamma, delta)) * pow(Type(2.0), Type(0.5) * delta - Type(1.0)) * distfun::mygammafn((delta + Type(1.0))/Type(2.0));
        return bc_mom;
    }

    template<class Type>
    Type std_aparch_moment(Type gamma, Type delta, Type shape){
        Type bc_mom = (pow(shape - Type(2.0),delta/Type(2.0)) * distfun::mygammafn((shape - delta)/Type(2.0)) *
            distfun::mygammafn((delta + Type(1.0))/Type(2.0)))/(distfun::mygammafn(shape/Type(2.0)) * Type(2.0) * sqrt(M_PI)) *
            (pow(Type(1.0) - gamma, delta) + pow(Type(1.0) + gamma,delta));
        return bc_mom;
    }

    template<class Type>
    Type ged_aparch_moment(Type gamma, Type delta, Type shape){
        Type lam = pow(Type(2.0), Type(-2.0)/shape) * distfun::mygammafn(Type(1.0)/shape) * (Type(1.0)/distfun::mygammafn(Type(3.0)/shape));
        Type lampower = sqrt(lam);
        Type bc_mom = (pow(Type(1.0) - gamma, delta) + pow(Type(1.0) + gamma, delta)) * ((distfun::mygammafn((delta + Type(1.0))/shape) *
            pow(lampower, delta) * pow(Type(2.0), (delta/shape) - Type(1.0)))/distfun::mygammafn(Type(1.0)/shape));
        return bc_mom;
    }

    // box-cox moment
    template<class Type>
    Type aparch_moment_func(Type gamma, Type delta, Type skew, Type shape, Type lambda, int dclass) {
        Type out = 0.0;
        switch(dclass){
        case 1:
            out = norm_aparch_moment(gamma, delta);
            break;
        case 2:
            out = std_aparch_moment(gamma, delta, shape);
            break;
        case 3:
            out = snorm_aparch_moment(gamma, delta, skew);
            break;
        case 4:
            out = sstd_aparch_moment(gamma, delta, skew, shape);
            break;
        case 5:
            out = ged_aparch_moment(gamma, delta, shape);
            break;
        case 6:
            out = sged_aparch_moment(gamma, delta, skew, shape);
            break;
        case 7:
            out = nig_aparch_moment(gamma, delta, skew, shape);
            break;
        case 8:
            out = gh_aparch_moment(gamma, delta, skew, shape, lambda);
            break;
        case 9:
            out = jsu_aparch_moment(gamma, delta, skew, shape);
            break;
        case 10:
            out = ghst_aparch_moment(gamma, delta, skew, shape);
            break;
        default:
            out = norm_aparch_moment(gamma, delta);
        }
        return out;
    }
}
