##-------------------------------------------------------##
## Auxiliary functions to fit censored longitudinal data ##
##-------------------------------------------------------##

## Normal multivariate density for censored data ##
dmvNCens <-function(cc, ni, y, tt,  X, U, D, beta, phi1, phi2, sigma2){
  N <- length(ni) # num. of subjects;
  n <- sum(ni) # total num. of observations;
  p <- dim(X)[2] # num. of beta parameters;
  q <- dim(U)[2] # num. of alpha parameters;
  
  dens1 <- c()
  D1 <- as.matrix(D) 
  for (i in 1:N){
    mu1 <- matrix(0,nrow= ni[i],ncol = 1)
    Omega1 <- matrix(0,nrow= ni[i],ncol = ni[i])
    E1 <- matrix(0,nrow= ni[i],ncol = ni[i])
    Sigma1 <- matrix(0,nrow= ni[i],ncol = ni[i])
    
    cc1 <- cc[(sum(ni[1:i-1])+1) : (sum(ni[1:i]))]
    y1 <- y[(sum(ni[1:i-1])+1) : (sum(ni[1:i]))]
    X1 <- matrix(X[(sum(ni[1:i-1])+1) : (sum(ni[1:i])), ],ncol=p)
    U1 <- matrix(U[(sum(ni[1:i-1])+1) : (sum(ni[1:i])), ],ncol=q)
    tt1 <- tt[(sum(ni[1:i-1])+1) : (sum(ni[1:i]))]
    mu1 <- X1%*%beta
    Omega1 <- DecM(tt1,phi1,phi2,sigma2)
    E1 <- matrix(Omega1/sigma2,nrow = ni[i],ncol = ni[i])
    Sigma1 <- (Omega1+(U1)%*%D1%*%t(U1))
    Sigma1 <- (Sigma1 + t(Sigma1))/2
    
    if(sum(cc1)==0){
      dens1[i] <- dmvnorm(x= as.vector(y1),mean =  as.vector(mu1), sigma = Sigma1)
    }
    if(sum(cc1)>0){
      if(sum(cc1)==ni[i]){
        dens1[i]<- pmnorm(x=as.vector(y1),mean=as.vector(mu1),varcov=Sigma1)
      }else{
        muc <- mu1[cc1==1]+Sigma1[cc1==1,cc1==0]%*%solve(Sigma1[cc1==0,cc1==0])%*%
          (y1[cc1==0]-mu1[cc1==0])
        Sc <- Sigma1[cc1==1,cc1==1]-Sigma1[cc1==1,cc1==0]%*%
          solve(Sigma1[cc1==0,cc1==0])%*%Sigma1[cc1==0,cc1==1]
        Sc <- (Sc + t(Sc))/2
        dens1[i] <- dmnorm(as.vector(y1[cc1==0]),as.vector(mu1[cc1==0]),Sigma1[cc1==0,cc1==0])*(pmnorm(as.vector(y1[cc1==1]),as.vector(muc),Sc))
      }
    }
  }
  return(dens1)
}
#---------------------------------------------#

## Mixture normal density for censored data ##
mixdmvNCens <- function(cc, ni, y, tt, X, U, D, beta, phi1, phi2, sigma2, pii){
  N <- length(ni) # num. of subjects;
  n <- sum(ni) # total num. of observations;
  p <- dim(X)[2] # num. of beta parameters;
  q <- dim(U)[2] # num. of alpha parameters;
  g <- dim(beta)[2]
  
  dens2 <- c(rep(0,N))
  
  for (G in 1:g) {
    dens2 <- dens2 + pii[G]*dmvNCens(cc=cc, ni=ni, y=y, tt=tt, X=X, U=U, D=D[[G]], beta=beta[,G], phi1=phi1[G], phi2=phi2[G], sigma2=sigma2[G])
  }
  return(dens2)
}
#---------------------------------------------#

## Normal multivariate density for censored data with uncorrelated structure ##
dmvNCens_unc <-function(cc, ni, y, tt,  X, U, D, beta, sigma2){
  N <- length(ni) # num. of subjects;
  n <- sum(ni) # total num. of observations;
  p <- dim(X)[2] # num. of beta parameters;
  q <- dim(U)[2] # num. of alpha parameters;
  
  dens1 <- c(rep(0,N))
  D1 <- as.matrix(D) 
  
  for (i in 1:N){
    mu1 <- matrix(0,nrow= ni[i],ncol = 1)
    Omega1 <- matrix(0,nrow= ni[i],ncol = ni[i])
    E1 <- matrix(0,nrow= ni[i],ncol = ni[i])
    Sigma1 <- matrix(0,nrow= ni[i],ncol = ni[i])
    
    cc1 <- cc[(sum(ni[1:i-1])+1) : (sum(ni[1:i]))]
    y1 <- y[(sum(ni[1:i-1])+1) : (sum(ni[1:i]))]
    X1 <- matrix(X[(sum(ni[1:i-1])+1) : (sum(ni[1:i])), ],ncol=p)
    U1 <- matrix(U[(sum(ni[1:i-1])+1) : (sum(ni[1:i])), ],ncol=q)
    tt1 <- tt[(sum(ni[1:i-1])+1) : (sum(ni[1:i]))]
    mu1 <- X1%*%beta
    Omega1 <- sigma2*diag(ni[i])
    E1 <- diag(ni[i])
    Sigma1 <- (Omega1+(U1)%*%D1%*%t(U1))
    Sigma1 <- (Sigma1 + t(Sigma1))/2
    
    if(sum(cc1)==0){
      dens1[i] <- dmvnorm(x= as.vector(y1),mean =  as.vector(mu1), sigma = Sigma1)
    }
    if(sum(cc1)>0){
      if(sum(cc1)==ni[i]){
        dens1[i]<- pmnorm(x=as.vector(y1),mean=as.vector(mu1),varcov=Sigma1)
      }else{
        muc <- mu1[cc1==1]+Sigma1[cc1==1,cc1==0]%*%solve(Sigma1[cc1==0,cc1==0])%*%
          (y1[cc1==0]-mu1[cc1==0])
        Sc <- Sigma1[cc1==1,cc1==1]-Sigma1[cc1==1,cc1==0]%*%
          solve(Sigma1[cc1==0,cc1==0])%*%Sigma1[cc1==0,cc1==1]
        Sc <- (Sc + t(Sc))/2
        dens1[i] <- dmnorm(as.vector(y1[cc1==0]),as.vector(mu1[cc1==0]),Sigma1[cc1==0,cc1==0])*(pmnorm(as.vector(y1[cc1==1]),as.vector(muc),Sc))
      }
    }
  }
  return(dens1)
}
#---------------------------------------------#

## Mixture normal density for censored data with uncorrelated structure ##
mixdmvNCens_unc <- function(cc, ni, y, tt, X, U, D, beta, sigma2, pii){
  N <- length(ni) # num. of subjects;
  n <- sum(ni) # total num. of observations;
  p <- dim(X)[2] # num. of beta parameters;
  q <- dim(U)[2] # num. of alpha parameters;
  g <- dim(beta)[2]
  
  dens2 <- 0
  
  for (G in 1:g) {
    dens2 <- dens2 + pii[G]*dmvNCens_unc(cc=cc, ni=ni, y=y, tt=tt, X=X, U=U, D=D[[G]], beta=beta[,G], sigma2=sigma2[G])
  }
  return(dens2)
}
#---------------------------------------------#

## Truncated moment for normal ##
Mnormtr <- function(u=c(0,0),S=diag(2),qc=c(1,2)) {
  
  nic=length(u)
  
  if (nic==1) {
    qq <- (1/sqrt(S))*(-qc+u)
    R<-1
    alpha <- pnorm(-qq)
    dd <- dnorm(-qq)
    H <- qq*dd
    EX <- (1/alpha)*dd   # a vector with a length of nic
    EXX <- 1+1/alpha*H
    varX <- EXX-EX^2
    Eycens <- -sqrt(S)*EX+u
    varyic<- varX*S
    E2yy<-varyic+Eycens^2
  }
  else {
    qq <- diag(1/sqrt(diag(S)))%*%(-qc+u)
    R <-  diag(1/sqrt(diag(S)))%*%S%*%diag(1/sqrt(diag(S)))
    alpha <- pmvnorm(upper=as.vector(-qq), corr=R)
    #print(qq)
    dd <- rep(0, nic)   #derivative vector
    
    for (j in 1:nic){
      V <- R[-j, -j, drop=F]-R[-j,j, drop=F]%*%R[j,-j, drop=F]
      nu <- -qq[-j]+R[-j,j, drop=F]%*%qq[j]
      dd[j] <- dnorm(-qq[j])*pmvnorm(upper=as.vector(nu), sigma=V)
    }
    
    H <- matrix(rep(0, nic*nic), nrow=nic)
    RH <- matrix(rep(0, nic*nic), nrow=nic)
    
    if(nic==2){
      H[1,2] <- H[2,1] <- dmvnorm(as.vector(-qq[c(1, 2)]),sigma=matrix(c(1, R[1,2], R[2,1], 1), nrow=2))
      #sigma==R since qq is standardized
      RH[1,2] <- RH[2,1] <- R[1,2]*H[1,2]
    }
    else {
      for( s in 1:(nic-1)){
        for (t in (s+1):nic){
          invR <- solve(R[c(s,t), c(s,t), drop=F])
          nu <- -qq[-c(s,t)]+R[-c(s,t), c(s,t), drop=F]%*%invR%*%qq[c(s,t),,drop=F]
          V <-  R[-c(s,t), -c(s,t), drop=F]- R[-c(s,t), c(s,t), drop=F]%*%invR%*%R[c(s,t), -c(s,t), drop=F]
          H[s,t] <- H[t,s] <- pmvnorm(upper=as.vector(nu), sigma=V)*dmvnorm(as.vector(-qq[c(s, t)]),sigma=matrix(c(1, R[s,t], R[t,s], 1), nrow=2))
          RH[s,t] <- RH[t,s] <- R[s,t]*H[s,t]
        }
      }
    }
    
    h <- qq*dd-apply(RH, 1, sum)
    diag(H) <- h
    EX <- (1/alpha)*R%*%dd   # a vector with a length of nic
    EXX <- R+1/alpha*R%*%H%*%R
    varX <- EXX-EX%*%t(EX)
    Eycens <- -diag(sqrt(diag(S)))%*%EX+u
    varyic <- diag(sqrt(diag(S)))%*%varX%*%diag(sqrt(diag(S)))
    E2yy <- varyic+Eycens%*%t(Eycens)
  }
  
  return(list(Ey=Eycens,Eyy=E2yy,Vary=varyic))
  
}
#---------------------------------------------#

## Correlation Matrix ##
# DEC Matrix
DecM <-function(tt,phi1,phi2,sigma2){
  aux_DecM <- (abs(outer(tt, tt, "-")))^phi2
  diag(aux_DecM) <-0
  Dec <- sigma2*(phi1^aux_DecM)
  return(Dec)
}
#---------------------------------------------#

# DEC matrix Derivative
DevEiAr1<-function(tt,phi1,phi2,sigma2){
  # rho = phi1    e gamma = phi2
  if(phi2<=0.0000001)
  {
    r <- length(tt)
    devR_phi1 <- matrix(1,nrow=r,ncol=r)
    diag(devR_phi1) <- 0
    devR_phi2 <- matrix(0,nrow=r,ncol=r)
  }
  else
  {
    func1 <- function(x,y){((abs(x-y))^phi2)*phi1^((abs(x-y))^phi2-1)}
    H1 <- (outer(tt, tt, func1))
    diag(H1) <- 0
    func2 <- function(x,y){((abs(x-y))^phi2)*log(abs(x-y))*log(phi1)*phi1^((abs(x-y))^phi2)}
    H2 <- (outer(tt, tt, func2))
    diag(H2) <- 0
    devR_phi1 <- H1
    devR_phi2 <- H2
  }
  
  obj.out <- list(devR_phi1 = devR_phi1, devR_phi2 = devR_phi2)
  return(obj.out)
}
#---------------------------------------------#

## Derivatives of the matrix of random effects with respect to each component
DerD<-function(M){
  
  m1<-dim(M)[1]
  m2<-dim(M)[2]  
  d<-list()
  for(h in 1:m1){
    d[[h]]<-list()
    for(k in 1:(m2+1-h)){
      d[[h]][[k]]<-matrix(0,m1,m2)
      if(k==1){d[[h]][[k]][h,h]<-1}   
      else{
        d[[h]][[k]][h,h+(k-1)]<-d[[h]][[k]][h+(k-1),h]<-1}
    }
  }
  
  return(d=d)
  
}
#---------------------------------------------#

## Estimate phi1 and phi2 for DEC ##
FCi1 <- function(phiG,beta,sigma2,tt,ubi,ubbi,uybi,uyyi,uyi,X,U,ni,Zij)
{
  phi11 <- phiG[1]  
  phi22 <- phiG[2]
  N <- length(ni)
  n<-sum(ni)
  p<-dim(X)[2]
  q<-dim(U)[2]
  aux1<- as.vector(Zij)
  
  soma <- c(0)
  
  for (i in 1:N){
    ub<-ubi[(((i-1)*q)+1) : (i*q), i]
    ubb<-ubbi[(((i-1)*q)+1) : (i*q), (((i-1)*q)+1) : (i*q)]
    uyb<-uybi[(sum(ni[1:i-1])+1) : (sum(ni[1:i])),(((i-1)*q)+1) : (i*q)]
    uyy<-uyyi[(sum(ni[1:i-1])+1) : (sum(ni[1:i])),(sum(ni[1:i-1])+1) : (sum(ni[1:i]))]
    uy<-uyi[(sum(ni[1:i-1])+1) : (sum(ni[1:i])),i]
    X1 <- matrix(X[(sum(ni[1:i-1])+1) : (sum(ni[1:i])),  ],ncol=p)
    U1 <- matrix(U[(sum(ni[1:i-1])+1) : (sum(ni[1:i])) ,  ],ncol=q)
    tt1 <- tt[(sum(ni[1:i-1])+1) : (sum(ni[1:i])) ]
    mup <- X1%*%beta  
    Cii <- DecM(tt1,phi11,phi22,sigma2)
    soma<- soma - 0.5*aux1[i]*log(det(Cii))-0.5*(sum(aux1[i]*diag(uyy%*%solve(Cii)))-(aux1[i]*t(uy)%*%solve(Cii)%*%mup)-(aux1[i]*t(mup)%*%solve(Cii)%*%uy)-sum(aux1[i]*diag(solve(Cii)%*%((uyb)%*%t(U1))))-sum(aux1[i]*diag(solve(Cii)%*%((uyb)%*%t(U1))))
                                                 +(aux1[i]*t(mup)%*%solve(Cii)%*%U1%*%ub)+(aux1[i]*t(ub)%*%t(U1)%*%solve(Cii)%*%mup)+(aux1[i]*t(mup)%*%solve(Cii)%*%mup)+sum(aux1[i]*diag(ubb%*%t(U1)%*%solve(Cii)%*%U1)))
  }
  return(as.vector(-soma))
}
#---------------------------------------------#

## phi2 fixed for AR  ##
FCi_phi1 <-function(phiG,phi22,beta,sigma2,tt,ubi,ubbi,uybi,uyyi,uyi,X,U,ni,Zij)
{
  phi11<-phiG
  N <- length(ni)
  n<-sum(ni)
  p<-dim(X)[2]
  q<-dim(U)[2]
  
  aux1<- as.vector(Zij)
  soma <- c(0)
  
  for (i in 1:N){
    ub<-ubi[(((i-1)*q)+1) : (i*q), i]
    ubb<-ubbi[(((i-1)*q)+1) : (i*q), (((i-1)*q)+1) : (i*q)]
    uyb<-uybi[(sum(ni[1:i-1])+1) : (sum(ni[1:i])),(((i-1)*q)+1) : (i*q)]
    uyy<-uyyi[(sum(ni[1:i-1])+1) : (sum(ni[1:i])),(sum(ni[1:i-1])+1) : (sum(ni[1:i]))]
    uy<-uyi[(sum(ni[1:i-1])+1) : (sum(ni[1:i])),i]
    X1 <- matrix(X[(sum(ni[1:i-1])+1) : (sum(ni[1:i])),  ],ncol=p)
    U1 <- matrix(U[(sum(ni[1:i-1])+1) : (sum(ni[1:i])) ,  ],ncol=q)
    tt1 <- tt[(sum(ni[1:i-1])+1) : (sum(ni[1:i])) ]
    mup <- X1%*%beta  
    Cii <- DecM(tt1,phi11,phi22,sigma2)
    
    soma<- soma - 0.5*aux1[i]*log(det(Cii))-0.5*(sum(diag(aux1[i]*uyy%*%solve(Cii)))-aux1[i]*t(uy)%*%solve(Cii)%*%mup-aux1[i]*t(mup)%*%solve(Cii)%*%uy-sum(diag(aux1[i]*solve(Cii)%*%((uyb)%*%t(U1))))-sum(diag(aux1[i]*solve(Cii)%*%((uyb)%*%t(U1))))
                                                 +aux1[i]*t(mup)%*%solve(Cii)%*%U1%*%ub+aux1[i]*t(ub)%*%t(U1)%*%solve(Cii)%*%mup+aux1[i]*t(mup)%*%solve(Cii)%*%mup+sum(diag(aux1[i]*ubb%*%t(U1)%*%solve(Cii)%*%U1)))
    
  }
  
  return(-soma)
}
#---------------------------------------------#