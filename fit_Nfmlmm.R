
source("structure_fit_functions.R")

Cens.fmlmEM <- function(cc, y, ni, tt, type, X ,U , g, type_structure, error = 0.00001, iter.max = 300){ 
  
  N <- length(ni) # num. of subjects;
  n <- sum(ni) # total num. of observations;
  p <- dim(X)[2] # num. of beta parameters;
  q <- dim(U)[2] # num. of alpha parameters;
  
  start.time <- Sys.time()  
  
  if(length(g) == 0) stop("g is not specified correctly.\n")
  
  # Initail values
  
  # To estimate Zij and pii    
  k.iter.max <- 50
  n.start    <- 1
  algorithm  <- "Hartigan-Wong"
  
  # begin initial values 
  if(g==1){pii <- rep(1/g,g)}
  if(g > 1){
    init <- kmeans(y,g,k.iter.max,n.start,algorithm)
    pii  <- as.vector(init$size/sum(ni)) 
  }
  
  mark1 <- lm(y~-1+X)
  
  beta <- matrix(rep(mark1$coefficients,g),ncol = g)
  
  mark2 <- summary(mark1)
  
  sigma2 <- rep((mark2$sigma)^2,g)

  D <- sD <- list() 
  
  for(j in 1:g){
    D[[j]] <- 0.1*diag(q)
    sD[[j]] <- solve(D[[j]])
  }
  
  structure_correlation = list("dec" = "dec", "ar" = "Phi1", 
                               "sym" = "Phi1", "ma" = "Phi1", "unc" = "unc")
  
  Struc = structure_correlation[[type_structure]]
  
  funclist[[Struc]](cc, y, ni, tt, X ,U , g, pii, type_structure, 
                    N, n, p, q, beta, sigma2, D, sD, start.time)
  
}# function
