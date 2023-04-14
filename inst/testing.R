library(xts)
library(TMB)
data(sp500ret, package = "rugarch")
data(dmbp, package = "rugarch")

devtools::clean_dll("src/")

pkgbuild::compile_dll()
devtools::document()
devtools::load_all()

y <- as.xts(sp500ret)
N <- nrow(dmbp)
yall <- as.xts(dmbp[,1], as.Date(1:N, origin = "1970-01-01"))
y <- as.xts(dmbp[1:(N - 100),1], as.Date(1:(N - 100), origin = "1970-01-01"))
vreg <- as.xts(dmbp[1:(N - 100),2], as.Date(1:(N - 100), origin = "1970-01-01"))
vreg_new <- as.xts(dmbp[((N - 100) + 1):N,2], as.Date(((N - 100) + 1):N, origin = "1970-01-01"))

spec <- garch_modelspec(y, constant = TRUE, model = "garch", init = "unconditional", vreg = vreg,
                        multiplicative = FALSE,
                        backcast_lambda = 0.9, sample_n = 10, variance_targeting = F, distribution = "jsu")
spec$parmatrix[parameter == "omega", lower := -0.2]
object <- estimate(spec, solver = "nloptr", control = nloptr_fast_options(maxeval = 100, trace = T))
sqrt(unconditional(object))
object$variance_target_summary

p <- predict(object, h = 100, newvreg = vreg_new, nsim = 2000)



plot(as.zoo(rbind(tail(sigma(object),1500), p$sigma)))
lines(as.zoo(p$sigma), col = 2)
abline(h = sqrt(unconditional(object)), col = 3)
lines(as.zoo(abs(tail(yall,100))), col = 4)

plot(abs(y), type = "l")
lines(sigma(object), col = 2)
summary(object)

print(summary(object, include_persistence = FALSE), format = "flextable")

halflife(object)

mod <- rugarch::ugarchfit(rugarch::ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = TRUE),
                          variance.model = list(variance.targeting = FALSE,
                                                model = "sGARCH", external.regressors = vreg), distribution.model = "jsu"),y, solver = "gosolnp", fit.control = list(stationarity = TRUE))
rugarch::persistence(mod)

cfx <- rugarch::coef(mod)
names(cfx)[5] <- "kappa1"
cbind(coef(object), cfx[names(coef(object))])
c(logLik(object), rugarch::likelihood(mod))
c(logLik(object) > rugarch::likelihood(mod))



benchmark <- fcp_benchmark()
lre_test(coef(object), benchmark[,"coefficient"])
lre_test(sqrt(diag(vcov(object))), benchmark[,"std_error_h"])
lre_test(sqrt(diag(vcov(object, type = "OP"))), benchmark[,"std_error_op"])
lre_test(sqrt(diag(vcov(object, type = "QMLE"))), benchmark[,"std_error_qmle"])

object2 <- tsfilter(object, y[1971:1971])

tail(sigma(object),2)
tail(sigma(object2),3)

spec <- garch_modelspec(y, constant = TRUE, variance_targeting = F, distribution = "norm")
object <- estimate(spec, solver = "nloptr", control = nloptr_default_options(trace = T))


coef(object)


spec <- garch_modelspec(y, variance_targeting = F, distribution = "norm", vreg = xts(dmbp[,2], index(y)), multiplicative = TRUE)
spec$parmatrix[parameter == "omega", lower := -10]
spec$parmatrix[parameter == "vreg1", lower := -20]

object <- estimate(spec, solver = "nloptr", control = nloptr_default_options(trace = T))
coef(object)

spec <- garch_modelspec(y, variance_targeting = F, distribution = "norm", vreg = xts(dmbp[,2], index(y)), multiplicative = FALSE)
object <- estimate(spec, solver = "nloptr", control = nloptr_default_options(trace = T))

sqrt(diag(vcov(object, type = "H")))
sqrt(diag(vcov(object, type = "OPG")))
sqrt(diag(vcov(object, type = "QMLE")))
sqrt(diag(vcov(object, type = "NW")))


control = list(trace = 1)

library(rugarch)
xspec <- ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = T), variance.model = list(garchOrder = c(1,1), variance.targeting = FALSE),
           distribution.model = "norm")
setfixed(xspec) <- as.list(coef(object))
mod<-ugarchfit(xspec,y[1:1970], solver = "solnp")

p <- ugarchforecast(xspec, data = y[1:1970], n.ahead = 20, n.roll = 0)


as.numeric(p@forecast$sigmaFor)

debugonce(rugarch:::.nsgarchforecast)

p <- ugarchforecast(xspec, data = y[1:1970], n.ahead = 20)
cbind(coef(object), coef(mod))

as.numeric(p@forecast$sigmaFor)

mod<-ugarchfit(ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = FALSE),
                          variance.model = list(variance.targeting = FALSE), distribution.model = "jsu"),y, solver = "solnp")


mod<-ugarchfit(ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = FALSE),
                          variance.model = list(variance.targeting = FALSE), distribution.model = "norm"),y, solver = "solnp")


cbind(sol$solution, coef(mod))
likelihood(mod)
L$fun(coef(mod), env)
L$tmb$fn(sol$par)

sol <- nlminb(start = cc, objective = tsgarch_env$fun, gradient = tsgarch_env$grad,
              lower = spec$parmatrix[estimate == 1]$lower, upper = spec$parmatrix[estimate == 1]$upper, control = list(trace = 1),
              tsgarch_env = tsgarch_env)

tsgarch_env$fun(spec$parmatrix[estimate == 1]$value, tsgarch_env)

tsgarch_env$tmb$report()$z


scores <- numDeriv::jacobian(func = score_function, x = scaled_solution$solution$par, env = scaled_solution$env)
H <- scaled_solution$env$tmb$he()
G <- diag(as.numeric(scaled_solution$env$tmb$gr()))
n <- NROW(scores)
V <- 1/n * (G %*% H %*% G)

jacobian(garch_ineq, x = optimal_pars, env = env)


f <- function(pars)
{
    ff <- function(x) djsu(x, mu = 0, sigma = 1, skew = pars[1], shape = pars[2])
    pars[3] + pars[4] * integrate(ff, -Inf, 0)$value
}

jacobian(f, x = c(0, 2, 0.3, 0.5))

jacobian(f, c(0,0,0))


eqs1 <- s$equation[[1]]
df <- data.frame(formula = eqs1)

cf <- s$coefficients
#empty <- data.table(term = as.character("."), Estimate = as.numeric(NULL), `Std. Error` = as.numeric(NULL), `t value` =  as.numeric(NULL),`Pr(>|t|)` = as.numeric(NULL))
print(flextable(cf) |> theme_apa() |> set_caption(caption = "GARCH Model Estimation") |>
    add_footer_lines(values = paste0("LogLik : ", round(-object$loglik,2))) |>
    add_footer_lines(values = paste0("AIC : ", round(AIC(object), 2))) |>
    add_footer_lines(values = " ") |>
    add_footer_lines(values = "Equation : ") |> append_chunks(i = 4, j = 1, part = "footer", as_equation(paste0("\\newline",s$equation[[1]],"\\newline",s$equation[[2]]))) |>
    compose(i = 1, j = 1, as_paragraph(as_chunk(''))) |> append_chunks(i=1,j=1, as_equation("\\mu")) |>
    compose(i = 2, j = 1, as_paragraph(as_chunk(''))) |> append_chunks(i=2,j=1, as_equation("\\omega")) |>
    compose(i = 3, j = 1, as_paragraph(as_chunk(''))) |> append_chunks(i=3,j=1, as_equation("\\alpha_1")) |>
    compose(i = 4, j = 1, as_paragraph(as_chunk(''))) |> append_chunks(i=4,j=1, as_equation("\\beta_1")), preview = "pdf")


# test an ets model and pass residuals to garch and predict and compare distribution
library(tsets)
library(tsdatasets)
y <- tsdatasets::electricload[,2]
ets_spec <- ets_modelspec(y, model = "AAA", frequency = 12)
ets_mod <- estimate(ets_spec, solver = "nlminb", autodiff = TRUE)
summary(ets_mod)
r <- residuals(ets_mod, raw = T)
spec <- garch_modelspec(r, constant = FALSE, variance_targeting = F, distribution = "norm")
#spec$parmatrix[parameter == "alpha1", lower := 1e-6]
control <- nloptr_default_options(trace = T)
control$maxeval <- 2000
object <- estimate(spec, solver = "nloptr", control = control)
summary(object)
plot(sigma(object))


mod <- rugarch::ugarchfit(rugarch::ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = FALSE), variance.model = list(garchOrder = c(1,1))),r)


########################################
# egarch
spec <- garch_modelspec(y[1:1970], constant = TRUE, model = "egarch", init = "unconditional",
                        backcast_lambda = 0.9, sample_n = 10, variance_targeting = T, distribution = "jsu")
control <- nloptr_fast_options(trace = T)
control$maxeval <- 500
#object <- spec
object <- estimate(spec, solver = "nloptr", control = control)
x <- summary(object)
print(summary(object))

f <- tsfilter(object, y = y[1971:1972])

f <- tsfilter(f, y = y[1973:1973])

plot(sigma(object) - sigma(f)[1:1970])

tail(sigma(f))

tail(sigma(object))

library(rugarch)
mod <- rugarch::ugarchfit(rugarch::ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = TRUE), variance.model = list(variance.targeting = T, model = "eGARCH", garchOrder = c(1,1)),
                                                    distribution.model = "ghst"),y)
cbind(coef(object), rugarch::coef(mod)[names(coef(object))])
c(-object$loglik, mod@fit$LLH)
c(-object$loglik > mod@fit$LLH)

plot(head(sigma(object), 100))
lines(head(rugarch::sigma(mod), 100), col = 2)

