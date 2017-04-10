library(waitup)
curve(from=0, to=10, expr=waitup::gamma_pdf(x=x, alpha=2, beta=1), col='green')
abline(v=qgamma(p=.5, shape=2, scale=1), col='green')
curve(from=0, to=10, expr=waitup:::generalized_gamma_pdf_3(x=x, alpha=2, beta=1, nu=.8), col='red', add=TRUE)
abline(v=flexsurv::qgengamma.orig(p=.5, shape=2, scale=1, k=2/.8), col='red')


