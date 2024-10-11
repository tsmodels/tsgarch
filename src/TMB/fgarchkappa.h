namespace fgarchkappa {
    // jsu
    template<class Float>
    struct struct_fgarch_jsu {
        typedef Float Scalar; // Required by integrate
        Float gamma, eta, delta, skew, shape;    // Data
        Float operator() (Float u) {
            Float ans = 0;
            const int give_log = 0;
            ans += pow(fabs(u - eta) - gamma * (u - eta), delta) * distfun::fwd_jsu_std(u, skew, shape, give_log);
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
    Float eval_fgarch_jsu(Float gamma, Float eta, Float delta, Float skew, Float shape) {
        struct_fgarch_jsu<Float> f = {gamma, eta, delta, skew, shape};
        return f.marginal();
    }

    TMB_BIND_ATOMIC(fgarchjsu, 11111, eval_fgarch_jsu(x[0], x[1], x[2], x[3], x[4]))

    template<class Type>
    Type jsu_fgarch_moment(Type gamma, Type eta, Type delta, Type skew, Type shape) {
            vector<Type> args(6); // Last index reserved for derivative order
            args << gamma, eta, delta, skew, shape, 0;
            return fgarchkappa::fgarchjsu(CppAD::vector<Type>(args))[0];
    }

    //sged
    template<class Float>
    struct struct_fgarch_sged {
        typedef Float Scalar; // Required by integrate
        Float gamma, eta, delta, skew, shape;    // Data
        Float operator() (Float u) {
            Float ans = 0;
            const int give_log = 0;
            ans += pow(fabs(u - eta) - gamma * (u - eta), delta) * distfun::fwd_sged_std(u, skew, shape, give_log);
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
    Float eval_fgarch_sged(Float gamma, Float eta, Float delta, Float skew, Float shape) {
        struct_fgarch_sged<Float> f = {gamma, eta, delta, skew, shape};
        return f.marginal();
    }

    TMB_BIND_ATOMIC(fgarchsged, 11111, eval_fgarch_sged(x[0], x[1], x[2], x[3], x[4]))

    template<class Type>
    Type sged_fgarch_moment(Type gamma, Type eta, Type delta, Type skew, Type shape) {
        vector<Type> args(6); // Last index reserved for derivative order
        args << gamma, eta, delta, skew, shape, 0;
        return fgarchkappa::fgarchsged(CppAD::vector<Type>(args))[0];
    }

    //sstd
    template<class Float>
    struct struct_fgarch_sstd {
        typedef Float Scalar; // Required by integrate
        Float gamma, eta, delta, skew, shape;    // Data
        Float operator() (Float u) {
            Float ans = 0;
            const int give_log = 0;
            ans += pow(fabs(u - eta) - gamma * (u - eta), delta) * distfun::fwd_sstd_std(u, skew, shape, give_log);
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
    Float eval_fgarch_sstd(Float gamma, Float eta, Float delta, Float skew, Float shape) {
        struct_fgarch_sstd<Float> f = {gamma, eta, delta, skew, shape};
        return f.marginal();
    }

    TMB_BIND_ATOMIC(fgarchsstd, 11111, eval_fgarch_sstd(x[0], x[1], x[2], x[3], x[4]))

    template<class Type>
    Type sstd_fgarch_moment(Type gamma, Type eta, Type delta, Type skew, Type shape) {
        vector<Type> args(6); // Last index reserved for derivative order
        args << gamma, eta, delta, skew, shape, 0;
        return fgarchkappa::fgarchsstd(CppAD::vector<Type>(args))[0];
    }

    //snorm
    template<class Float>
    struct struct_fgarch_snorm {
        typedef Float Scalar; // Required by integrate
        Float gamma, eta, delta, skew;    // Data
        Float operator() (Float u) {
            Float ans = 0;
            const int give_log = 0;
            ans += pow(fabs(u - eta) - gamma * (u - eta), delta) * distfun::fwd_snorm_std(u, skew, give_log);
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
    Float eval_fgarch_snorm(Float gamma, Float eta, Float delta, Float skew) {
        struct_fgarch_snorm<Float> f = {gamma, eta, delta, skew};
        return f.marginal();
    }
    TMB_BIND_ATOMIC(fgarchsnorm, 1111, eval_fgarch_snorm(x[0], x[1], x[2], x[3]))
    template<class Type>
    Type snorm_fgarch_moment(Type gamma, Type eta, Type delta, Type skew) {
        vector<Type> args(5); // Last index reserved for derivative order
        args << gamma, eta, delta, skew, 0;
        return fgarchkappa::fgarchsnorm(CppAD::vector<Type>(args))[0];
    }

    //nig
    template<class Float>
    struct struct_fgarch_nig {
        typedef Float Scalar; // Required by integrate
        Float gamma, eta, delta, skew, shape;    // Data
        Float operator() (Float u) {
            Float ans = 0;
            const int give_log = 0;
            ans += pow(fabs(u - eta) - gamma * (u - eta), delta) * distfun::fwd_nig_std(u, skew, shape, give_log);
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
    Float eval_fgarch_nig(Float gamma, Float eta, Float delta, Float skew, Float shape) {
        struct_fgarch_nig<Float> f = {gamma, eta, delta, skew, shape};
        return f.marginal();
    }

    TMB_BIND_ATOMIC(fgarchnig, 11111, eval_fgarch_nig(x[0], x[1], x[2], x[3], x[4]))

    template<class Type>
    Type nig_fgarch_moment(Type gamma, Type eta, Type delta, Type skew, Type shape) {
        vector<Type> args(6); // Last index reserved for derivative order
        args << gamma, eta, delta, skew, shape, 0;
        return fgarchkappa::fgarchnig(CppAD::vector<Type>(args))[0];
    }


      template<class Float>
      struct struct_fgarch_ghst {
          typedef Float Scalar; // Required by integrate
          Float gamma, eta, delta, skew, shape;    // Data
          Float operator() (Float u) {
              Float ans = 0;
              const int give_log = 0;
              ans += pow(fabs(u - eta) - gamma * (u - eta), delta) * distfun::fwd_ghst_std(u, skew, shape, give_log);
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
      Float eval_fgarch_ghst(Float gamma, Float eta, Float delta, Float skew, Float shape) {
          struct_fgarch_ghst<Float> f = {gamma, eta, delta, skew, shape};
          return f.marginal();
      }

      TMB_BIND_ATOMIC(fgarchghst, 11111, eval_fgarch_ghst(x[0], x[1], x[2], x[3], x[4]))

      template<class Type>
      Type ghst_fgarch_moment(Type gamma, Type eta, Type delta, Type skew, Type shape) {
          vector<Type> args(6); // Last index reserved for derivative order
          args << gamma, eta, delta, skew, shape, 0;
          return fgarchkappa::fgarchghst(CppAD::vector<Type>(args))[0];
      }

      template<class Float>
      struct struct_fgarch_gh {
          typedef Float Scalar; // Required by integrate
          Float gamma, eta, delta, skew, shape, lambda;    // Data
          Float operator() (Float u) {
              Float ans = 0;
              const int give_log = 0;
              ans += pow(fabs(u - eta) - gamma * (u - eta), delta) * distfun::fwd_gh_std(u, skew, shape, lambda, give_log);
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
      Float eval_fgarch_gh(Float gamma, Float eta, Float delta, Float skew, Float shape, Float lambda) {
          struct_fgarch_gh<Float> f = {gamma, eta, delta, skew, shape, lambda};
          return f.marginal();
      }

      TMB_BIND_ATOMIC(fgarchgh, 111111, eval_fgarch_gh(x[0], x[1], x[2], x[3], x[4], x[5]))

      template<class Type>
      Type gh_fgarch_moment(Type gamma, Type eta, Type delta, Type skew, Type shape, Type lambda) {
          vector<Type> args(7); // Last index reserved for derivative order
          args << gamma, eta, delta, skew, shape, lambda, 0;
          return fgarchkappa::fgarchgh(CppAD::vector<Type>(args))[0];
      }

    // we don't have closed form solutions for the non skewed distributions

    //norm
    template<class Float>
    struct struct_fgarch_norm {
        typedef Float Scalar; // Required by integrate
        Float gamma, eta, delta;    // Data
        Float operator() (Float u) {
            Float ans = 0;
            const int give_log = 0;
            ans += pow(fabs(u - eta) - gamma * (u - eta), delta) * distfun::norm_std(u, give_log);
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
    Float eval_fgarch_norm(Float gamma, Float eta, Float delta) {
        struct_fgarch_norm<Float> f = {gamma, eta, delta};
        return f.marginal();
    }

    TMB_BIND_ATOMIC(fgarchnorm, 111, eval_fgarch_norm(x[0], x[1], x[2]))

    template<class Type>
    Type norm_fgarch_moment(Type gamma, Type eta, Type delta) {
        vector<Type> args(4); // Last index reserved for derivative order
        args << gamma, eta, delta, 0;
        return fgarchkappa::fgarchnorm(CppAD::vector<Type>(args))[0];
    }

    //std
    template<class Float>
    struct struct_fgarch_std {
        typedef Float Scalar; // Required by integrate
        Float gamma, eta, delta, shape;    // Data
        Float operator() (Float u) {
            Float ans = 0;
            const int give_log = 0;
            ans += pow(fabs(u - eta) - gamma * (u - eta), delta) * distfun::fwd_student_std(u, shape, give_log);
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
    Float eval_fgarch_std(Float gamma, Float eta, Float delta, Float shape) {
        struct_fgarch_std<Float> f = {gamma, eta, delta, shape};
        return f.marginal();
    }

    TMB_BIND_ATOMIC(fgarchstd, 1111, eval_fgarch_std(x[0], x[1], x[2], x[3]))

    template<class Type>
    Type std_fgarch_moment(Type gamma, Type eta, Type delta, Type shape) {
        vector<Type> args(5); // Last index reserved for derivative order
        args << gamma, eta, delta, shape, 0;
        return fgarchkappa::fgarchstd(CppAD::vector<Type>(args))[0];
    }

    //ged
    template<class Float>
    struct struct_fgarch_ged {
        typedef Float Scalar; // Required by integrate
        Float gamma, eta, delta, shape;    // Data
        Float operator() (Float u) {
            Float ans = 0;
            const int give_log = 0;
            ans += pow(fabs(u - eta) - gamma * (u - eta), delta) * distfun::ged_std(u, shape, give_log);
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
    Float eval_fgarch_ged(Float gamma, Float eta, Float delta, Float shape) {
        struct_fgarch_ged<Float> f = {gamma, eta, delta, shape};
        return f.marginal();
    }

    TMB_BIND_ATOMIC(fgarchged, 1111, eval_fgarch_ged(x[0], x[1], x[2], x[3]))

    template<class Type>
    Type ged_fgarch_moment(Type gamma, Type eta, Type delta, Type shape) {
        vector<Type> args(5); // Last index reserved for derivative order
        args << gamma, eta, delta, shape, 0;
        return fgarchkappa::fgarchged(CppAD::vector<Type>(args))[0];
    }


    // box-cox moment
    template<class Type>
    Type fgarch_moment_func(Type gamma, Type eta, Type delta, Type skew, Type shape, Type lambda, int dclass) {
        Type out = 0.0;
        switch(dclass){
        case 1:
            out = norm_fgarch_moment(gamma, eta, delta);
            break;
        case 2:
            out = std_fgarch_moment(gamma, eta, delta, shape);
            break;
        case 3:
            out = snorm_fgarch_moment(gamma, eta, delta, skew);
            break;
        case 4:
            out = sstd_fgarch_moment(gamma, eta, delta, skew, shape);
            break;
        case 5:
            out = ged_fgarch_moment(gamma, eta, delta, shape);
            break;
        case 6:
            out = sged_fgarch_moment(gamma, eta, delta, skew, shape);
            break;
        case 7:
            out = nig_fgarch_moment(gamma, eta, delta, skew, shape);
            break;
        case 8:
            out = gh_fgarch_moment(gamma, eta, delta, skew, shape, lambda);
            break;
        case 9:
            out = jsu_fgarch_moment(gamma, eta, delta, skew, shape);
            break;
        case 10:
            out = ghst_fgarch_moment(gamma, eta, delta, skew, shape);
            break;
        default:
            out = norm_fgarch_moment(gamma, eta, delta);
        }
        return out;
    }
}
