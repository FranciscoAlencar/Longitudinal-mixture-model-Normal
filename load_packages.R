##-------------------------------------------------------##
##      Packages to fit censored longitudinal data       ##
##-------------------------------------------------------##

packageurl <- "http://cran.r-project.org/src/contrib/Archive/mvtnorm/mvtnorm_1.0-8.tar.gz"
install.packages(packageurl, repos=NULL, type="source")

library("mvtnorm")

library("mnormt")

library("lmec")

library("tmvtnorm")

library("MASS")

library("optimx")

library("mixsmsn")
