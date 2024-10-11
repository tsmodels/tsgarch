likelihood_fun <- function(x, distribution = "norm", mu = 0, sigma = 1, skew = 0, shape = 0, lambda = 0)
{
    out <- sum(log(ddist(distribution = distribution, x = (x - mu)/sigma, mu = 0, sigma = 1, skew = skew, shape = shape, lambda = lambda)/sigma))
    return(out)
}
