code <- readLines('genes8PT_analytical_figures.R')
eval(parse(text=code[seq_len(grep('^# Verify closed-form',code)-1)]))
for (B in c(5000,25000)) {
 s <- .01; Ne <- 10000; x0 <- 1/(2*Ne)
 lc <- -log1p(-recombination_loss(B));lu <- -B*log1p(-mu)
 b <- log1p(s)-lc;a <- b-lu;xstar <- a/b
 times <- half_sweep(s,Ne,B,TRUE)
 x <- function(t) xstar*plogis(qlogis(x0/xstar)+a*t)
 stopifnot(abs(x(times$T)-.5)<1e-13)
 A_num <- integrate(function(t) 1-x(t),0,times$T,rel.tol=1e-10)$value
 T_num <- integrate(function(y) 1/(y*(a-b*y)),x0,.5,rel.tol=1e-10)$value
 stopifnot(abs(A_num-times$A)<1e-8,abs(T_num-times$T)<1e-8)
}
# Genetic boundary coordinates must be independent of bp-to-Morgan calibration.
reference <- axis_value(retention_bound(c(.01,.1,.25),1000,.95),'cM')
rho <- 3e-8
stopifnot(max(abs(reference-axis_value(retention_bound(c(.01,.1,.25),1000,.95),'cM')))<1e-12)
mu <- 0
for(B in c(1,1000,20000)) {
 stopifnot(abs(log_retention(.1,1000,B,TRUE)-log_retention(.1,1000,B,FALSE))<1e-12)
}
cat('Independent mutation quadrature, zero-mutation limit and cM invariance passed.\n')