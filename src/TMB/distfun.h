namespace distfun{
    // fwd_function_name : fwd are the functions rewritten for tapeless forward mode

    // add custom besselk function to allow scaling (exponential scaling flag = 2)
    TMB_BIND_ATOMIC(bessel_k2, 11, atomic::bessel_utils::bessel_k(x[0], x[1], 2.))

    template<class Type>
    Type scaled_besselK(Type x, Type nu)
    {
        Type ans;
        CppAD::vector<Type> tx(3);
        tx[0] = x;
        tx[1] = nu;
        tx[2] = 0;
        ans = bessel_k2(tx)[0];
        return ans;
    }
    VECTORIZE2_tt(scaled_besselK)

    template<class Type>
    Type fwd_besselK(Type x, Type nu)
    {
        Type ans = atomic::bessel_utils::bessel_k(x, nu, 1.);
        return ans;
    }
    VECTORIZE2_tt(fwd_besselK)

    template<class Type>
    Type fwd_scaled_besselK(Type x, Type nu)
    {
        Type ans = atomic::bessel_utils::bessel_k(x, nu, 2.);
        return ans;
    }
    VECTORIZE2_tt(fwd_scaled_besselK)


    template <class Type>
    Type mygammafn(Type x)
    {
        Type out = exp(lgamma(x));
        return out;
    }
    VECTORIZE1_t(mygammafn)

    template <class Type>
    Type norm_std(Type x, int give_log)
    {
        Type pdf;
        pdf = dnorm(x, Type(0.0), Type(1.0), give_log);
        return pdf;
    }
    VECTORIZE1_t(norm_std)

    template <class Type>
    Type student_std(Type x, Type shape, int give_log)
    {
        Type pdf;
        if (shape <= Type(2.0)) {
            pdf = 1E12;
        } else {
            Type s = sqrt(shape/(shape - Type(2.0)));
            pdf = dt(x * s, shape, 0) * s;
        }
        if (give_log == 1) pdf = log(pdf);
        return pdf;
    }
    VECTORIZE3_tti(student_std)


    template <class Type>
    Type jsu_std(Type x, Type skew, Type shape, int give_log)
    {
        Type rtau = Type(1.0)/shape;
        Type w = CppAD::CondExpLt(rtau, Type(0.0000001), Type(1.0), exp(rtau*rtau));
        Type omega = -skew * rtau;
        Type c = sqrt(Type(1.0)/(Type(0.5)*(w - Type(1.0))*(w*cosh(Type(2.0)*omega)+Type(1.0))));
        Type z = (x - (c * sqrt(w) * sinh(omega)))/c;
        Type r = -skew + log(z + sqrt(z*z+Type(1.0)))/rtau;
        Type pdf = -log(c) - log(rtau) - Type(0.5) * log( z * z + Type(1.0)) - Type(0.5) * log(Type(2.0) * M_PI) - Type(0.5) * r * r;
        if(give_log == 0) {
            pdf = exp(pdf);
        }
        return pdf;
    }
    VECTORIZE4_ttti(jsu_std)

    template <class Type>
    Type sstd_std(Type x, Type skew, Type shape, int give_log)
    {
        Type a = Type(1.0)/Type(2.0);
        Type b = shape/Type(2.0);
        Type beta = exp((lgamma(a) - lgamma(a+b)) + lgamma(b));
        Type m1 = Type(2.0) * sqrt(shape - Type(2.0))/(shape - Type(1.0))/beta;
        Type mu = m1 * (skew - Type(1.0)/skew);
        Type sigma = sqrt((Type(1.0) - m1 * m1) * (skew * skew + Type(1.0)/(skew * skew)) + Type(2.0) * m1 * m1 - Type(1.0));
        Type z = x * sigma + mu;
        Type xxi_tmp = CppAD::CondExpLt(z, Type(0.0), Type(1.0)/skew, skew);
        Type xxi = CppAD::CondExpEq(z, Type(0.0), Type(1.0), xxi_tmp);
        Type g = Type(2.0)/(skew+Type(1.0)/skew);
        Type pdf = g * student_std(z/xxi,shape,0) * sigma;
        if(give_log == 1) pdf = log(pdf);
        return pdf;
    }
    VECTORIZE4_ttti(sstd_std)

    template <class Type>
    Type ged_std(Type x, Type shape, int give_log)
    {
        Type lambda = sqrt(pow(Type(1.0)/Type(2.0),Type(2.0)/shape)*mygammafn(Type(1.0)/shape)/mygammafn(Type(3.0)/shape));
        Type g = shape/(lambda*(pow(Type(2.0),Type(1.0)+(Type(1.0)/shape)))*mygammafn(Type(1.0)/shape));
        Type pdf = g * exp(Type(-0.5)*pow(fabs(x/lambda),shape));
        if(give_log == 1) pdf = log(pdf);
        return pdf;
    }
    VECTORIZE3_tti(ged_std)

    template <class Type>
    Type sged_std(Type x, Type skew, Type shape, int give_log)
    {
        Type lambda = sqrt(pow(Type(1.0)/Type(2.0),(Type(2)/shape))*mygammafn(Type(1.0)/shape)/mygammafn(Type(3.0)/shape));
        Type m1 = pow(Type(2.0), (Type(1.0)/shape))*lambda*mygammafn(Type(2.0)/shape)/mygammafn(Type(1.0)/shape);
        Type mu = m1*(skew-Type(1.0)/skew);
        Type skew2 = skew * skew;
        Type m12 = m1 * m1;
        Type sigma = sqrt((Type(1.0) - m12)*(skew2 + Type(1.0)/skew2) + Type(2.0) * m12 - Type(1.0));
        Type z = x*sigma + mu;
        Type xxi_tmp = CppAD::CondExpLt(z, Type(0.0), Type(1.0)/skew, skew);
        Type xxi = CppAD::CondExpEq(z, Type(0.0), Type(1.0), xxi_tmp);
        Type g = Type(2.0)/(skew + Type(1.0)/skew);
        Type pdf = g * ged_std(z/xxi, shape, 0) * sigma;
        if(give_log == 1) pdf = log(pdf);
        return pdf;
    }
    VECTORIZE4_ttti(sged_std)

    template <class Type>
    Type snorm_std(Type x, Type skew, int give_log)
    {
        Type m1 = Type(2.0)/sqrt(Type(2.0)*M_PI);
        Type m12 = m1 * m1;
        Type xi2 = skew*skew;
        Type mu = m1*(skew-Type(1.0)/skew);
        Type sigma = sqrt((Type(1.0)-m12) * (xi2 + Type(1.0)/xi2) + Type(2.0) * m12 - Type(1.0));
        Type z = x*sigma+mu;
        Type xxi = CppAD::CondExpLt(z, Type(0.0), Type(1.0)/skew, skew);
        Type g = Type(2.0)/(skew + Type(1.0)/skew);
        Type pdf = g * dnorm(z/xxi,Type(0), Type(1.0), 0) * sigma;
        if(give_log == 1) pdf = log(pdf);
        return pdf;
    }
    VECTORIZE3_tti(snorm_std)


    template <class Type>
    Type kappagh(Type x, Type lambda)
    {
        Type kappa = CppAD::CondExpEq(lambda, Type(-0.5), Type(1.0)/x, ((scaled_besselK(x, lambda + Type(1.0)) / (scaled_besselK(x, lambda))) / x));
        return kappa;
    }
    VECTORIZE2_tt(kappagh)

    template <class Type>
    Type fwd_kappagh(Type x, Type lambda)
    {
        Type kappa = Type(0.0);
        Type tmp = lambda + Type(1.0);
        kappa = (fwd_scaled_besselK(x, tmp) /fwd_scaled_besselK(x, lambda)) / x;
        return kappa;
    }
    VECTORIZE2_tt(fwd_kappagh)

    template <class Type>
    Type deltakappagh(Type x, Type lambda)
    {
        Type dkappa = kappagh(x, Type(lambda + Type(1.0))) - kappagh(x, lambda);
        return dkappa;
    }
    VECTORIZE2_tt(deltakappagh)

    template <class Type>
    Type fwd_deltakappagh(Type x, Type lambda)
    {
        Type tmp = lambda + Type(1.0);
        Type dkappa = fwd_kappagh(x, tmp) - fwd_kappagh(x, lambda);
        return dkappa;
    }
    VECTORIZE2_tt(fwd_deltakappagh)

    template <class Type>
    Type nig_std(Type x, Type skew, Type shape, int give_log)
    {
        Type rho2 = Type(1.0) - skew * skew;
        Type zeta2 = shape * shape;
        Type alpha = zeta2 * kappagh(shape, Type(-0.5))/rho2;
        alpha = alpha * ( Type(1.0) + skew * skew * zeta2 * deltakappagh(shape, Type(-0.5))/rho2);
        alpha = sqrt(alpha);
        Type beta = alpha * skew;
        Type delta = shape / ( alpha * sqrt(rho2) );
        Type d = delta*delta;
        Type mu = -beta * (d) * kappagh(shape, Type(-0.5));
        Type xm = x - mu;
        Type pdf =  -log(M_PI) + log(alpha) + log(delta) + log(besselK(alpha * sqrt(d + xm * xm), Type(1.0))) +
            delta*sqrt(alpha*alpha-beta*beta)+beta*xm-Type(0.5)*log(d+xm*xm);
        if(give_log == 0){
            pdf = exp(pdf);
        }
        return pdf;
    }
    VECTORIZE4_ttti(nig_std)

    template <class Type>
    Type fwd_nig_std(Type x, Type skew, Type shape, int give_log)
    {
        Type rho2 = Type(1.0) - skew * skew;
        Type zeta2 = shape * shape;
        Type alpha = zeta2 * fwd_kappagh(shape, Type(-0.5))/rho2;
        alpha = alpha * (Type(1.0) + skew * skew * zeta2 * fwd_deltakappagh(shape, Type(-0.5))/rho2);
        alpha = sqrt(alpha);
        Type beta = alpha * skew;
        Type delta = shape / ( alpha * sqrt(rho2) );
        Type d = delta*delta;
        Type mu = -beta * (d) * fwd_kappagh(shape, Type(-0.5));
        Type xm = x - mu;
        Type tmp = alpha * sqrt(d + xm * xm);
        Type pdf =  -log(M_PI) + log(alpha) + log(delta) + log(fwd_besselK(tmp, Type(1.0))) +
            delta*sqrt(alpha*alpha-beta*beta)+beta*xm-Type(0.5)*log(d+xm*xm);
        if(give_log == 0){
            pdf = exp(pdf);
        }
        return pdf;
    }
    VECTORIZE4_ttti(fwd_nig_std)

    template <class Type>
    Type ghst_std(Type x, Type skew, Type shape, int give_log)
    {
        // zero undefined
        Type shapeplus = CppAD::CondExpLt(fabs(shape), Type(1e-12), Type(1e-12), shape);
        Type shapem2 = shapeplus - Type(2.0);
        Type delta = sqrt(Type(1.0)/(((Type(2.0) * skew*skew)/(shapem2*shapem2*(shapeplus-Type(4.0)))) + (Type(1.0)/shapem2)));
        Type beta = skew/delta;
        Type beta2 = beta * beta;
        Type delta2 = delta * delta;
        Type mu = -( (beta * (delta*delta))/shapem2);
        Type xmu = x - mu;
        Type xmu2 = xmu * xmu;
        Type arg = sqrt(beta2 * (delta2 + xmu2));
        Type arg2 = (shapeplus + Type(1.0))/Type(2.0);
        Type pdf = ((Type(1.0) - shapeplus)/Type(2.0)) * log(Type(2.0)) + shapeplus * log(delta) +
            arg2 * log(fabs(beta)) + log(scaled_besselK(arg, arg2)) - sqrt(beta2 * (delta2 + xmu2)) +
            beta * xmu - lgamma(shapeplus/Type(2.0)) - log(M_PI)/Type(2.0) -
            arg2 * log(delta2 + xmu2)/Type(2.0);
        if(give_log == 0){
            pdf = exp(pdf);
        }
        return pdf;
    }
    VECTORIZE4_ttti(ghst_std)

    template <class Type>
    Type fwd_ghst_std(Type x, Type skew, Type shape, int give_log)
    {
        // zero undefined
        Type shapeplus = shape;
        if (fabs(shape)<Type(1e-12)) {
            shapeplus = Type(1e-12);
        }
        Type shapem2 = shapeplus - Type(2.0);
        Type delta = sqrt(Type(1.0)/( ((Type(2.0) * skew*skew)/(shapem2*shapem2*(shapeplus-Type(4.0)))) + (Type(1.0)/shapem2) ));
        Type beta = skew/delta;
        Type beta2 = beta * beta;
        Type delta2 = delta * delta;
        Type mu = -((beta * (delta*delta))/shapem2);
        Type xmu = x - mu;
        Type xmu2 = xmu * xmu;
        Type arg = sqrt(beta2 * (delta2 + xmu2));
        Type arg2 = (shapeplus + Type(1.0))/Type(2.0);
        Type tmp = fwd_scaled_besselK(arg, arg2);
        Type pdf = ((Type(1.0) - shapeplus)/Type(2.0)) * log(Type(2.0)) +
            shapeplus * log(delta) + arg2 * log(fabs(beta)) +
            log(tmp) - sqrt(beta2 * (delta2 + xmu2)) +
            beta * xmu - lgamma(shapeplus/Type(2.0)) -
            log(M_PI)/Type(2.0) - arg2 * log(delta2 + xmu2)/Type(2.0);
        if(give_log == 0){
            pdf = exp(pdf);
        }
        return pdf;
    }
    VECTORIZE4_ttti(fwd_ghst_std)

    template <class Type>
    Type gh(Type x, Type alpha, Type beta, Type delta, Type mu, Type lambda)
    {
        Type pdf = 0.0;
        if(alpha <= 0){
            return pdf;
        }
        if(delta <= 0){
            return pdf;
        }
        if(fabs(beta) >= alpha){
            return pdf;
        }
        Type alpha2 = alpha * alpha;
        Type beta2 = beta * beta;
        Type delta2 = delta * delta;
        Type arg = delta*sqrt(alpha2 - beta2);
        Type xmu = x - mu;
        Type xmu2 = xmu * xmu;
        Type a = (lambda/Type(2.0)) * log(alpha2 - beta2) - (log(sqrt(Type(2.0) * M_PI)) +
            (lambda - Type(0.5)) * log(alpha) + lambda * log(delta) +
            log(scaled_besselK(arg, lambda)) - arg);
        Type f = ((lambda - Type(0.5))/Type(2.0)) * log(delta2 + xmu2);
        arg = alpha * sqrt(delta2 + xmu2);
        Type k = log(scaled_besselK(arg, lambda - Type(0.5))) - arg;
        Type e = beta * xmu;
        pdf = exp(a + f + k + e);
        return pdf;
    }

    template <class Type>
    Type fwd_gh(Type x, Type alpha, Type beta, Type delta, Type mu, Type lambda)
    {
        Type pdf = 0.0;
        if(alpha <= 0){
            return pdf;
        }
        if(delta <= 0){
            return pdf;
        }
        if(fabs(beta) >= alpha){
            return pdf;
        }
        Type alpha2 = alpha * alpha;
        Type beta2 = beta * beta;
        Type delta2 = delta * delta;
        Type arg = delta*sqrt(alpha2 - beta2);
        Type xmu = x - mu;
        Type xmu2 = xmu * xmu;
        Type tmp1 = fwd_scaled_besselK(arg, lambda);
        Type a = (lambda/Type(2.0)) * log(alpha2 - beta2) - (log(sqrt(Type(2.0) * M_PI)) +
            (lambda - Type(0.5)) * log(alpha) + lambda * log(delta) + log(tmp1) - arg);
        Type f = ((lambda - Type(0.5))/Type(2.0)) * log(delta2 + xmu2);
        arg = alpha * sqrt(delta2 + xmu2);
        Type tmp2 = lambda - Type(0.5);
        Type tmp3 = fwd_scaled_besselK(arg, tmp2);
        Type k = log(tmp3) - arg;
        Type e = beta * xmu;
        pdf = exp(a + f + k + e);
        return pdf;
    }


    template <class Type>
    Type gh_std(Type x, Type skew, Type shape, Type lambda, int give_log)
    {
        Type rho2 = Type(1.0) - skew * skew;
        Type zeta2 = shape * shape;
        Type alpha_tmp = zeta2 * kappagh(shape, lambda)/rho2;
        Type alpha = sqrt(alpha_tmp * (Type(1.0) + skew * skew * zeta2 * deltakappagh(shape, lambda)/rho2));
        Type beta = alpha * skew;
        Type delta = shape / ( alpha * sqrt(rho2) );
        Type d = delta*delta;
        Type mu = -beta * (d) * kappagh(shape, lambda);
        Type pdf = gh(x, alpha, beta, delta, mu, lambda);
        if(give_log == 1){
            pdf = log(pdf);
        }
        return pdf;
    }
    VECTORIZE5_tttti(gh_std)

    template <class Type>
    Type fwd_gh_std(Type x, Type skew, Type shape, Type lambda, int give_log)
    {
        Type rho2 = Type(1.0) - skew * skew;
        Type zeta2 = shape * shape;
        Type alpha_tmp = zeta2 * fwd_kappagh(shape, lambda)/rho2;
        Type alpha = sqrt(alpha_tmp * (Type(1.0) + skew * skew * zeta2 * fwd_deltakappagh(shape, lambda)/rho2));
        Type beta = alpha * skew;
        Type delta = shape / ( alpha * sqrt(rho2) );
        Type d = delta*delta;
        Type mu = -beta * (d) * fwd_kappagh(shape, lambda);
        Type pdf = fwd_gh(x, alpha, beta, delta, mu, lambda);
        if(give_log == 1){
            pdf = log(pdf);
        }
        return pdf;
    }
    VECTORIZE5_tttti(fwd_gh_std)

    template<class Type>
    Type signbranch(Type x, Type u){
        Type out = u;
        if (x < Type(0.0)) {
            out = Type(1.0)/u;
        }
        if (x == Type(0.0)) {
            out = Type(1.0);
        }
        return out;
    }
    VECTORIZE1_t(signbranch)

    template <class Type>
    Type fwd_jsu_std(Type x, Type skew, Type shape, int give_log)
    {
        Type rtau = Type(1.0)/shape;
        Type w = exp(rtau * rtau);
        Type omega = -skew * rtau;
        Type c = sqrt(Type(1.0)/(Type(0.5)*(w - Type(1.0))*(w*cosh(Type(2.0)*omega)+Type(1.0))));
        Type z = (x - (c * sqrt(w) * sinh(omega)))/c;
        Type r = -skew + log(z + sqrt(z*z+Type(1.0)))/rtau;
        Type pdf = -log(c) - log(rtau) - Type(0.5) * log( z * z + Type(1.0)) - Type(0.5) * log(Type(2.0) * M_PI) - Type(0.5) * r * r;
        if(give_log == 0) {
            pdf = exp(pdf);
        }
        return pdf;
    }
    VECTORIZE4_ttti(fwd_jsu_std)

    template <class Type>
    Type fwd_student_std(Type x, Type shape, int give_log)
    {
        Type pdf;
        Type s = sqrt(shape/(shape - Type(2.0)));
        Type tmp = s*x;
        pdf = dt(tmp, shape, 0) * s;
        if (give_log == 1) pdf = log(pdf);
        return pdf;
    }
    VECTORIZE3_tti(fwd_student_std)

    template <class Type>
    Type fwd_sstd_std(Type x, Type skew, Type shape, int give_log)
    {
        Type a = Type(1.0)/Type(2.0);
        Type b = shape/Type(2.0);
        Type ab = a + b;
        Type tmp1 = lgamma(Type(1.0) * a);
        Type tmp2 = lgamma(Type(1.0) * ab);
        Type tmp3 = lgamma(Type(1.0) * b);
        Type tmp4 = (tmp1 - tmp2) + tmp3;
        Type beta = exp(tmp4);
        Type m1 = Type(2.0) * sqrt(shape - Type(2.0))/(shape - Type(1.0))/beta;
        Type mu = m1 * (skew - Type(1.0)/skew);
        Type sigma = sqrt((Type(1.0) - m1 * m1) * (skew * skew + Type(1.0)/(skew * skew)) + Type(2.0) * m1 * m1 - Type(1.0));
        Type z = x * sigma + mu;
        Type xxi = signbranch(z, skew);
        Type g = Type(2.0)/(skew+Type(1.0)/skew);
        Type tmp = z/xxi;
        Type pdf = g * fwd_student_std(tmp,shape,0) * sigma;
        if(give_log == 1) pdf = log(pdf);
        return pdf;
    }
    VECTORIZE4_ttti(fwd_sstd_std)


    template <class Type>
    Type fwd_sged_std(Type x, Type skew, Type shape, int give_log)
    {
        Type lambda = sqrt(pow(Type(1.0)/Type(2.0),(Type(2)/shape))*mygammafn(Type(1.0)/shape)/mygammafn(Type(3.0)/shape));
        Type m1 = pow(Type(2.0), (Type(1.0)/shape))*lambda*mygammafn(Type(2.0)/shape)/mygammafn(Type(1.0)/shape);
        Type mu = m1*(skew-Type(1.0)/skew);
        Type skew2 = skew * skew;
        Type m12 = m1 * m1;
        Type sigma = sqrt((Type(1.0) - m12)*(skew2 + Type(1.0)/skew2) + Type(2.0) * m12 - Type(1.0));
        Type z = x*sigma + mu;
        Type xxi = signbranch(z, skew);
        Type g = Type(2.0)/(skew + Type(1.0)/skew);
        Type tmp = z/xxi;
        Type pdf = g * ged_std(tmp, shape, 0) * sigma;
        if(give_log == 1) pdf = log(pdf);
        return pdf;
    }
    VECTORIZE4_ttti(fwd_sged_std)

    template <class Type>
    Type fwd_snorm_std(Type x, Type skew, int give_log)
    {
        Type m1 = Type(2.0)/sqrt(Type(2.0)*M_PI);
        Type m12 = m1 * m1;
        Type xi2 = skew*skew;
        Type mu = m1*(skew-Type(1.0)/skew);
        Type sigma = sqrt((Type(1.0)-m12) * (xi2 + Type(1.0)/xi2) + Type(2.0) * m12 - Type(1.0));
        Type z = x*sigma+mu;
        Type xxi = signbranch(z, skew);
        Type g = Type(2.0)/(skew + Type(1.0)/skew);
        Type tmp = z/xxi;
        Type pdf = g * dnorm(tmp,Type(0), Type(1.0), 0) * sigma;
        if(give_log == 1) pdf = log(pdf);
        return pdf;
    }
    VECTORIZE3_tti(fwd_snorm_std)

    template <class Type>
    Type distlike(Type x, Type skew, Type shape, Type lambda, int dclass)
    {
        Type out = 0.0;
        const int give_log = 0;
        switch(dclass){
        case 1:
            out = norm_std(x,give_log);
            break;
        case 2:
            out = student_std(x, shape, give_log);
            break;
        case 3:
            out = snorm_std(x, skew, give_log);
            break;
        case 4:
            out = sstd_std(x, skew, shape, give_log);
            break;
        case 5:
            out = ged_std(x, shape, give_log);
            break;
        case 6:
            out = sged_std(x, skew, shape, give_log);
            break;
        case 7:
            out = nig_std(x, skew, shape, give_log);
            break;
        case 8:
            out = gh_std(x, skew, shape, lambda, give_log);
            break;
        case 9:
            out = jsu_std(x, skew, shape, give_log);
            break;
        case 10:
            out = ghst_std(x, skew, shape, give_log);
            break;
        default:
            out = norm_std(x,give_log);
        }
        return(out);
    }
    VECTORIZE5_tttti(distlike)

}
