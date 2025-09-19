x=1
theta0=2
v=1
r=1


lambda2 <- function(t, v, theta0, xbar,r) {
  max(v*(theta0+v*t-xbar),0)+r
}
lambda <- Vectorize(lambda2, "t")

print(integrate(lambda, 0, 1, v = v, theta0 = theta0, xbar = x, r=r))

ftn2 <- function(t, v, theta0, xbar, r) {
  lambda(t, v, theta0, xbar, r)*exp(-integrate(lambda, 0, t, v = v, theta0 = theta0, xbar = xbar, r=r)$value)
}
ft <- Vectorize(ftn2, "t")

Ft2 <- function(t, v, theta0, xbar, r) {
  if(t == 0)
    return(0)
  integrate(ft, 0, t, v = v, theta0 = theta0, xbar = xbar, r=r)$value
}
Ft <- Vectorize(Ft2, "t")

Ft.inv2 <- function(z, v, theta0, xbar, r) {
  uniroot(\(t) z-Ft(t, v, theta0, xbar, r), 0, 100)
}
Ft.inv <- Vectorize(Ft.inv2, "z")

# curve(lambda(x, v = v, theta0 = theta0, xbar = xbar, r=r), from = -2, to = 2)
# curve(ft(x, v = v, theta0 = theta0, xbar = xbar, r=r), from = -2, to = 2)
# curve(Ft(x, v = v, theta0 = theta0, xbar = xbar, r=r), from = 0, to = 2)

# Likelihood ratio
curve(ft(x, v = v, theta0 = theta0, xbar = 1,r=r) /
      ft(x, v = v, theta0 = theta0, xbar = 0, r=r), from = 0, to = 5)

# tradeoff function for test H0: xbar1 vs H1: xbar2
# TODO: take account of v
beta2 <- function(alpha, xbar1, xbar2, v, theta0, r) {
  if(xbar1 < xbar2) {
    # f(|xbar1) / f(|xbar2)
    # => monotone decreasing
    t <- uniroot(\(t) (1-alpha)-Ft(t, v, theta0, xbar1, r),c(0, 100))$root
    Ft(t, v, theta0, xbar2, r)
  } else {
    # if xbar1 <= xbar2
    # f(|xbar1) / f(|xbar2)
    # => monotone increasing
    t <- uniroot(\(t) alpha-Ft(t, v, theta0, xbar1, r), c(0, 100))$root
    1-Ft(t, v, theta0, xbar2, r)
  }
}
beta <- Vectorize(beta2, "alpha")

library(fdp)

print(fdp(
  zigzag1 = fdp_point(beta(alpha, xbar1 = 1, xbar2 = 0, v = v, theta0 = theta0, r=1)),
  zigzag5 = fdp_point(beta(alpha, xbar1 = 1, xbar2 = 0, v = v, theta0 = theta0, r=5)),
  zigzag10 = fdp_point(beta(alpha, xbar1 = 1, xbar2 = 0, v = v, theta0 = theta0, r=10)),
  zigzag20 = fdp_point(beta(alpha, xbar1 = 1, xbar2 = 0, v = v, theta0 = theta0, r=20)),
  zigzag30 = fdp_point(beta(alpha, xbar1 = 1, xbar2 = 0, v = v, theta0 = theta0, r=30)),
  gdp(1)))

