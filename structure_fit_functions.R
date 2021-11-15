source("aux_functions.R")
dec_function = function(cc, y, ni, tt, X ,U , g, pii, type_structure, 
                        N, n, p, q, beta, sigma2, D, sD, start.time,
                        error = 0.00001, iter.max = 300){

  cat("Runnig model with", type_structure, "correlation structure ... \n")
  
  phi1 <- c(rep(0.5,g))
  
  phi2 <- c(rep(2,g))
  
  theta <- c(beta[,1],sigma2[1],D[[1]][upper.tri(D[[1]], diag = T)],phi1[1],phi2[1])
  
  n_alphas <- length(D[[1]][upper.tri(D[[1]], diag = T)])
  
  d <- g*length(theta) + (g-1)
  
  criteria <- 1
  count <- 0
  lkante <- 1
    
  while((criteria > error) && (count <= iter.max)){
    
    si.beta <- si.sigma <- si.phi1 <- si.phi2 <- si.pii <- si <- si.g <- NULL
    
    MI <- matrix(0,g*sum(c(p,1,1,1,1,q*(1+ q)/2)) - 1 ,g*sum(c(p,1,1,1,1,q*(1+ q)/2)) - 1)     
    
    count <- count + 1
    
    print(count)
    
    Zij <- matrix(0, N, g)  
    
    ubii <- list()
    
    ubi <- matrix(0,N*q,N)
    ubbi <- matrix(0,N*q,N*q)
    uybi <- matrix(0,n,N*q)
    uyyi <- matrix(0,n,n)
    uyi <- matrix(0,n,N)
    
    lik <- matrix(0,N,g)
    
    soma5 <- matrix(0,N,g)
    
    for(j in 1:g){
      soma1 <- matrix(0,q,q)
      soma2 <- c(rep(0,g)) # It is referred to the sigma2
      soma3 <- matrix(0,p,p)
      soma4 <- matrix(0,p,g) # It is referred to the second sum beta's
      
      if(sqrt(sum((order(pii)-g:1)*(order(pii)-g:1))) != 0 ){
        betai  <- beta
        sigma2i <- sigma2
        phi1i <- phi1
        phi2i <- phi2
        for(l in 1:g){
          beta[,l]    <- betai[,(order(pii, decreasing = TRUE)[l])]
          sigma2[l] <- sigma2i[(order(pii, decreasing = TRUE)[l])]
          phi1[l] <- phi1i[(order(pii, decreasing = TRUE)[l])]
          phi2[l] <- phi2i[(order(pii, decreasing = TRUE)[l])]
        }
        pii <- pii[order(pii, decreasing = TRUE)]
      }
      
      d1 <- dmvNCens(cc=cc, ni=ni, y=y, tt= tt, X=X, U=U, D=D[[j]], beta=beta[,j], phi1= phi1[j], phi2=phi2[j], sigma2=sigma2[j])
      if(length(which(d1 == 0)) > 0) d1[which(d1 == 0)] <- .Machine$double.xmin 
      
      d2 <- mixdmvNCens(cc=cc, ni=ni, y=y, tt= tt, X=X, U=U, D=D, beta=beta, phi1= phi1, phi2=phi2, sigma2=sigma2, pii=pii)      
      if(length(which(d2 == 0)) > 0) d2[which(d2 == 0)] <- .Machine$double.xmin 
      
      Zij[,j] <- (pii[j]*d1)/(d2)
      
      pii[j] <- (1/N)*sum(Zij[,j])
      
      for (i in 1:N){
        cci <- cc[(sum(ni[1 :i-1])+1) : (sum(ni[1:i]))]
        yi <- y[(sum(ni[1:i-1])+1) : (sum(ni[1:i]))]
        Xi <- matrix(X[(sum(ni[1:i-1])+1) : (sum(ni[1:i])), ],ncol=p)
        Ui <- matrix(U[(sum(ni[1:i-1])+1) : (sum(ni[1:i])), ],ncol=q)
        tti <- tt[(sum(ni[1:i-1])+1) : (sum(ni[1:i]))]
        mu <- Xi%*%beta[,j]
        Omega <- DecM(tti,phi1[j],phi2[j],sigma2[j])
        E <- Omega/sigma2[j]
        Sigma <- (Omega+(Ui)%*%D[[j]]%*%t(Ui))
        Sigma <- (Sigma + t(Sigma))/2
        delta <- solve(sD[[j]]+(t(Ui)%*%solve(E)%*%(Ui*(1/sigma2[j]))))
        delta <- (delta + t(delta))/2
        
        if(sum(cci)==0){
          uy <- matrix(yi,ni[i],1)
          uyy <- yi%*%t(yi)
          ub <- delta%*%(t(Ui)*(1/sigma2[j]))%*%solve(E)%*%(uy-mu)
          ubb <- delta+(delta%*%(t(Ui)*((1/sigma2[j])^2))%*%solve(E)%*%(uyy-uy%*%t(mu) - mu%*%t(uy)+mu%*%t(mu))%*%solve(E)%*%Ui%*%delta)
          uyb <- (uyy-uy%*%t(mu))%*%solve(E)%*%(Ui*(1/sigma2[j]))%*%delta
          lik[i,j] <- dmvnorm(as.vector(yi),as.vector(mu),Sigma)
        }
        
        if(sum(cci)>=1){
          if(sum(cci)==ni[i]){
            muc <- Xi%*%beta[,j]
            Sc <- Sigma
            aux<- Mnormtr(muc,Sc,yi)
            uy<-aux$Ey
            uyy<- aux$Eyy
            ub <- delta%*%(t(Ui)*(1/sigma2[j]))%*%solve(E)%*%(uy-mu)
            ubb <- delta+(delta%*%(t(Ui)*((1/sigma2[j])^2))%*%solve(E)%*%(uyy-uy%*%t(mu) - mu%*%t(uy)+mu%*%t(mu))%*%solve(E)%*%Ui%*%delta)
            uyb <- (uyy-uy%*%t(mu))%*%solve(E)%*%(Ui*(1/sigma2[j]))%*%delta
          }
          else {
            muc <- Xi[cci==1,]%*%beta[,j]+Sigma[cci==1,cci==0]%*%
              solve(Sigma[cci==0,cci==0])%*%(yi[cci==0]-Xi[cci==0,]%*%beta[,j])
            Sc <- Sigma[cci==1,cci==1]-Sigma[cci==1,cci==0]%*%solve(Sigma[cci==0,cci==0])%*%Sigma[cci==0,cci==1]
            Sc <- (Sc + t(Sc))/2
            aux <- mtmvnorm(as.vector(muc), Sc, lower = rep(-Inf, length = length(muc)), 
                            upper = as.vector(yi[cci==1]))
            uy <- matrix(yi,ni[i],1)
            uy[cci==1] <- aux$tmean
            uyy <- matrix(0,ni[i],ni[i])
            uyy[cci==1,cci==1] <- aux$tvar
            uyy <- uyy+uy%*%t(uy)
            ub <- delta%*%(t(Ui)*(1/sigma2[j]))%*%solve(E)%*%(uy-mu)
            ubb <- delta+(delta%*%(t(Ui)*((1/sigma2[j])^2))%*%solve(E)%*%(uyy-uy%*%t(mu) - mu%*%t(uy)+mu%*%t(mu))%*%solve(E)%*%Ui%*%delta)
            uyb <- (uyy-uy%*%t(mu))%*%solve(E)%*%(Ui*(1/sigma2[j]))%*%delta
          }
        }
        
        soma1 <- soma1 + (Zij[i,j]*ubb)/sum(Zij[,j])          
        
        soma2[j] <- soma2[j] + (sum(Zij[i,j]*diag(uyy%*%solve(E))) - (Zij[i,j]*t(uy)%*%solve(E)%*%mu) - (Zij[i,j]*t(mu)%*%solve(E)%*%uy) - sum(Zij[i,j]*diag(solve(E)%*%((uyb)%*%t(Ui))))-sum(Zij[i,j]*diag(solve(E)%*%((uyb)%*%t(Ui))))
                                +Zij[i,j]*t(mu)%*%solve(E)%*%Ui%*%ub+Zij[i,j]*t(ub)%*%t(Ui)%*%solve(E)%*%mu
                                +Zij[i,j]*t(mu)%*%solve(E)%*%mu+sum(Zij[i,j]*diag(ubb%*%t(Ui)%*%solve(E)%*%Ui)))/(sum(ni[i]*Zij[,j]))
        
        soma3 <- soma3 + Zij[i,j]*(t(Xi)%*%solve(E)%*%Xi)
        
        soma4[,j] <- soma4[,j] + Zij[i,j]*(t(Xi)%*%solve(E)%*%(uy-Ui%*%ub))
        
        soma5[i,j] <- pii[j]*lik[i,j]
        
        ubi[(((i-1)*q)+1) : (i*q), i] <- ub
        
        ubbi[(((i-1)*q)+1) : (i*q), (((i-1)*q) + 1) : (i*q)] <- ubb
        
        uybi[(sum(ni[1:i-1])+1) : (sum(ni[1:i])),(((i-1)*q)+1) : (i*q)] <- uyb
        
        uyyi[(sum(ni[1:i-1])+1) : (sum(ni[1:i])),(sum(ni[1:i-1])+1) : (sum(ni[1:i]))] <- uyy
        
        uyi[(sum(ni[1:i-1])+1) : (sum(ni[1:i])),i] <- uy
        
        # Information matrix
        
        # pi
        si.pi <- (Zij[i,j]/as.numeric(pii[j])) - (Zij[i,g]/as.numeric(pii[g]))
        
        # beta
        si.beta <- (1/sigma2[j])*((t(Xi)%*%solve(E)%*%(Zij[i,j]*uy-Zij[i,j]*Ui%*%ub)) - 
                                    (Zij[i,j]*t(Xi)%*%solve(E)%*%Xi%*%beta[,j]))
        #sigmae
        si.sigma <- -(1/2)*(((Zij[i,j]*ni[i])/sigma2[j])-(1/sigma2[j]^2)*((sum(diag(Zij[i,j]*uyy%*%solve(E)))-Zij[i,j]*t(uy)%*%solve(E)%*%mu-Zij[i,j]*t(mu)%*%solve(E)%*%uy-sum(diag(Zij[i,j]*solve(E)%*%((uyb)%*%t(Ui))))-sum(diag(Zij[i,j]*solve(E)%*%((uyb)%*%t(Ui)))) +
                                                                             Zij[i,j]*t(mu)%*%solve(E)%*%Ui%*%ub+Zij[i,j]*t(ub)%*%t(Ui)%*%solve(E)%*%mu+Zij[i,j]*t(mu)%*%solve(E)%*%mu+sum(diag(Zij[i,j]*ubb%*%t(Ui)%*%solve(E)%*%Ui)))))
        
        # phi1 
        Dp <- DevEiAr1(tti,phi1[j],phi2[j],sigma2[j]) 
        Dpr_phi1 <- Dp$devR_phi1
        Dpr_phi2 <- Dp$devR_phi2
        dE1_phi1 <- sum(diag(solve(E)%*%Dpr_phi1))
        dE2_phi1 <- -(solve(E)%*%Dpr_phi1%*%solve(E))
        dE1_phi2 <- sum(diag(solve(E)%*%Dpr_phi2))
        dE2_phi2 <- -(solve(E)%*%Dpr_phi2%*%solve(E)) 
        si.phi1 <- - 0.5*Zij[i,j]*dE1_phi1 - (0.5/sigma2[j])*((sum(diag(Zij[i,j]*uyy%*%dE2_phi1))-Zij[i,j]*t(uy)%*%dE2_phi1%*%mu-Zij[i,j]*t(mu)%*%dE2_phi1%*%uy-sum(diag(Zij[i,j]*dE2_phi1%*%((uyb)%*%t(Ui))))-sum(diag(Zij[i,j]*dE2_phi1%*%((uyb)%*%t(Ui)))) +
                                                                 Zij[i,j]*t(mu)%*%dE2_phi1%*%Ui%*%ub+Zij[i,j]*t(ub)%*%t(Ui)%*%dE2_phi1%*%mu+Zij[i,j]*t(mu)%*%dE2_phi1%*%mu+sum(diag(Zij[i,j]*ubb%*%t(Ui)%*%dE2_phi1%*%Ui))))      
        si.phi2 <- - 0.5*Zij[i,j]*dE1_phi2 - (0.5/sigma2[j])*((sum(diag(Zij[i,j]*uyy%*%dE2_phi2))-Zij[i,j]*t(uy)%*%dE2_phi2%*%mu-Zij[i,j]*t(mu)%*%dE2_phi2%*%uy-sum(diag(Zij[i,j]*dE2_phi2%*%((uyb)%*%t(Ui))))-sum(diag(Zij[i,j]*dE2_phi2%*%((uyb)%*%t(Ui)))) +
                                                                 Zij[i,j]*t(mu)%*%dE2_phi2%*%Ui%*%ub+Zij[i,j]*t(ub)%*%t(Ui)%*%dE2_phi2%*%mu+Zij[i,j]*t(mu)%*%dE2_phi2%*%mu+sum(diag(Zij[i,j]*ubb%*%t(Ui)%*%dE2_phi2%*%Ui))))      
        
        # D
        D_der <- DerD(D[[j]])
        deralpha <- rep(0,n_alphas)
        md2<-dim(D[[j]])[1]  
        kont <- 0
        for(i1 in 1:md2){
          for(i2 in 1:(md2+1-i1)){
            kont <- kont+1
            di <- D_der[[i1]][[i2]]
            deralpha[kont] <- (-0.5)*sum(diag(Zij[i,j]*sD[[j]]%*%di-Zij[i,j]*sD[[j]]%*%di%*%sD[[j]]*ubb))     
          }
        }
        
        si <- rbind(si,c(si.beta,si.sigma,si.phi1,si.phi2,deralpha,si.pi))
        
      }
      
      si.g <- cbind(si.g,si)
      
      si <- NULL
      
      ubii[[j]] <- ubi
      
      beta[,j] <- solve(soma3)%*%(soma4[,j])
      sigma2[j] <- soma2[j]
      D[[j]] <- soma1
      sD[[j]] <- solve(D[[j]]) 
      
      phis <- optimx(c(phi1[j],phi2[j]), fn = FCi1, method = "nlminb",
                     beta = beta[,j],sigma2=sigma2[j],tt=tt,ubi=ubi,ubbi=ubbi,uybi=uybi,uyyi=uyyi,uyi=uyi,X=X,U=U,
                     ni=ni,Zij = Zij[,j], lower=c(0.01,0.01), upper=c(0.9,30),hessian=TRUE)
      phi1[j] <- phis$p1
      phi2[j] <- phis$p2
      
      pii[g] <- 1 - (sum(pii) - pii[g])
      
      
    }

    for(v in 1:N) MI <- MI + si.g[v,-g*sum(c(p,1,1,1,1,q*(1+q)/2))]%*%t(si.g[v,-g*sum(c(p,1,1,1,1,q*(1+ q)/2))])         
    
    lk <- sum(log(mixdmvNCens(cc=cc, ni=ni, y=y, tt= tt, X=X, U=U, D=D, beta=beta,phi1= phi1, phi2=phi2, sigma2=sigma2, pii=pii)))
    criteria <- abs((lk/lkante-1))
    lkante <- lk

  }# while
    
  end.time <- Sys.time()
  time.taken <- end.time - start.time
    
  aic <- -2*lk + 2*d
  bic <- -2*lk + log(N)*d
  edc <- -2*lk + 0.2*sqrt(N)*d
    
  out <- list(d=d,sigma2 = sigma2,pii = pii,phi1= phi1, phi2= phi2,beta = beta,D=D,ubi = ubii,Z = Zij,MI = MI,AIC=aic,BIC=bic,EDC=edc,loglik=lk,time = time.taken)
    
  return(out)     
    
}

Phi1_function = function(cc, y, ni, tt, X ,U , g, pii, type_structure, 
                         N, n, p, q, beta, sigma2, D, sD, start.time,
                         error = 0.00001, iter.max = 300){
  phi1 <- c(rep(0.5,g))
  
  Phi_2 = list("ar" = c(rep(1,g)), "sym" = c(rep(0,g)), "ma" = c(rep(Inf,g))) 
  
  cat("Runnig model with", type_structure, "correlation structure ... \n" )
  
  phi2 = Phi_2[[type_structure]]
  
  theta <- c(beta[,1],sigma2[1],D[[1]][upper.tri(D[[1]], diag = T)],phi1[1])
  n_alphas <- length(D[[1]][upper.tri(D[[1]], diag = T)])
  
  d <- g*length(theta) + (g-1)

  criteria <- 1
  count <- 0
  lkante <- 1
  
  while((criteria > error) && (count <= iter.max) ){
    
    si.beta <- si.sigma <- si.phi1 <- si.pii <- si <- si.g <- NULL
    
    MI <- matrix(0,g*sum(c(p,1,1,1,q*(1+ q)/2)) - 1 ,g*sum(c(p,1,1,1,q*(1+ q)/2)) - 1)   
    
    count <- count + 1
    
    print(count)
    
    Zij <- matrix(0, N, g)  
    
    ubii <- list()
    
    ubi <- matrix(0,N*q,N)
    ubbi <- matrix(0,N*q,N*q)
    uybi <- matrix(0,n,N*q)
    uyyi <- matrix(0,n,n)
    uyi <- matrix(0,n,N)
    
    lik <- matrix(0,N,g)
    
    soma5 <- matrix(0,N, g)
    
    
    for (j in 1:g){
    
      if( sqrt(sum((order(pii)-g:1)*(order(pii)-g:1))) != 0 ){
        betai  <- beta
        sigma2i <- sigma2
        phi1i <- phi1
        for(l in 1:g){
          beta[,l]    <- betai[,(order(pii, decreasing = TRUE)[l])]
          sigma2[l] <- sigma2i[(order(pii, decreasing = TRUE)[l])]
          phi1[l] <- phi1i[(order(pii, decreasing = TRUE)[l])]
        }
        pii <- pii[order(pii, decreasing = TRUE)]
      }  
      
      soma1 <- matrix(0,q,q)
      soma2 <- c(rep(0,g)) # It is referred to the sigma2
      soma3 <- matrix(0,p,p)
      soma4 <- matrix(0,p,g) # It is referred to the second sum beta's
      
      d1 <- dmvNCens(cc=cc, ni=ni, y=y, tt= tt, X=X, U=U, D=D[[j]], beta=beta[,j],phi1= phi1[j], phi2=phi2[j], sigma2=sigma2[j])
      if(length(which(d1 == 0)) > 0) d1[which(d1 == 0)] <- .Machine$double.xmin # de donde pertenece
      
      d2 <- mixdmvNCens(cc=cc, ni=ni, y=y, tt= tt, X=X, U=U, D=D, beta=beta,phi1= phi1, phi2=phi2, sigma2=sigma2, pii=pii)      
      if(length(which(d2 == 0)) > 0) d2[which(d2 == 0)] <- .Machine$double.xmin # de donde pertenece
      
      Zij[,j] <- (pii[j]*d1)/(d2)
      
      pii[j] <- (1/N)*sum(Zij[,j])
      
      for (i in 1:N){
    
        cci <- cc[(sum(ni[1 :i-1])+1) : (sum(ni[1:i]))]
        yi <- y[(sum(ni[1:i-1])+1) : (sum(ni[1:i]))]
        Xi <- matrix(X[(sum(ni[1:i-1])+1) : (sum(ni[1:i])), ],ncol=p)
        Ui <- matrix(U[(sum(ni[1:i-1])+1) : (sum(ni[1:i])), ],ncol=q)
        tti <- tt[(sum(ni[1:i-1])+1) : (sum(ni[1:i]))]
        mu <- Xi%*%beta[,j]
        Omega <- DecM(tti,phi1[j],phi2[j],sigma2[j])
        E <- Omega/sigma2[j]
        Sigma <- (Omega+(Ui)%*%D[[j]]%*%t(Ui))
        Sigma <- (Sigma + t(Sigma))/2
        delta <- solve(sD[[j]]+(t(Ui)%*%solve(E)%*%(Ui*(1/sigma2[j]))))
        delta <- (delta + t(delta))/2
        
        if(sum(cci)==0){
          uy <- matrix(yi,ni[i],1)
          uyy <- yi%*%t(yi)
          ub <- delta%*%(t(Ui)*(1/sigma2[j]))%*%solve(E)%*%(uy-mu)
          ubb <- delta+(delta%*%(t(Ui)*((1/sigma2[j])^2))%*%solve(E)%*%(uyy-uy%*%t(mu) - mu%*%t(uy)+mu%*%t(mu))%*%solve(E)%*%Ui%*%delta)
          uyb <- (uyy-uy%*%t(mu))%*%solve(E)%*%(Ui*(1/sigma2[j]))%*%delta
        }
        
        if(sum(cci)>=1){
          if(sum(cci)==ni[i]){
            muc <- Xi%*%beta[,j]
            Sc <- Sigma
            aux<- Mnormtr(muc,Sc,yi)
            uy<-aux$Ey
            uyy<- aux$Eyy
            ub <- delta%*%(t(Ui)*(1/sigma2[j]))%*%solve(E)%*%(uy-mu)
            ubb <- delta+(delta%*%(t(Ui)*((1/sigma2[j])^2))%*%solve(E)%*%(uyy-uy%*%t(mu) - mu%*%t(uy)+mu%*%t(mu))%*%solve(E)%*%Ui%*%delta)
            uyb <- (uyy-uy%*%t(mu))%*%solve(E)%*%(Ui*(1/sigma2[j]))%*%delta
          }
          else {
            muc <- Xi[cci==1,]%*%beta[,j]+Sigma[cci==1,cci==0]%*%
              solve(Sigma[cci==0,cci==0])%*%(yi[cci==0]-Xi[cci==0,]%*%beta[,j])
            Sc <- Sigma[cci==1,cci==1]-Sigma[cci==1,cci==0]%*%solve(Sigma[cci==0,cci==0])%*%Sigma[cci==0,cci==1]
            Sc <- (Sc + t(Sc))/2
            aux <- Mnormtr(muc,Sc,yi[cci==1])
            uy <- matrix(yi,ni[i],1)
            uy[cci==1] <- aux$Ey
            uyy <- matrix(0,ni[i],ni[i])
            uyy[cci==1,cci==1] <- aux$Vary
            uyy <- uyy+uy%*%t(uy)
            ub <- delta%*%(t(Ui)*(1/sigma2[j]))%*%solve(E)%*%(uy-mu)
            ubb <- delta+(delta%*%(t(Ui)*((1/sigma2[j])^2))%*%solve(E)%*%(uyy-uy%*%t(mu) - mu%*%t(uy)+mu%*%t(mu))%*%solve(E)%*%Ui%*%delta)
            uyb <- (uyy-uy%*%t(mu))%*%solve(E)%*%(Ui*(1/sigma2[j]))%*%delta
          }
        }
        
        soma1 <- soma1 + (Zij[i,j]*ubb)/sum(Zij[,j])          
        
        soma2[j] <- soma2[j] + (sum(Zij[i,j]*diag(uyy%*%solve(E))) - (Zij[i,j]*t(uy)%*%solve(E)%*%mu) - (Zij[i,j]*t(mu)%*%solve(E)%*%uy) - sum(Zij[i,j]*diag(solve(E)%*%((uyb)%*%t(Ui))))-sum(Zij[i,j]*diag(solve(E)%*%((uyb)%*%t(Ui))))
                                +Zij[i,j]*t(mu)%*%solve(E)%*%Ui%*%ub+Zij[i,j]*t(ub)%*%t(Ui)%*%solve(E)%*%mu
                                +Zij[i,j]*t(mu)%*%solve(E)%*%mu+sum(Zij[i,j]*diag(ubb%*%t(Ui)%*%solve(E)%*%Ui)))/(sum(ni[i]*Zij[,j]))
        
        soma3 <- soma3 + Zij[i,j]*(t(Xi)%*%solve(E)%*%Xi)
        
        soma4[,j] <- soma4[,j] + Zij[i,j]*(t(Xi)%*%solve(E)%*%(uy-Ui%*%ub))
        
        soma5[i,j] <- pii[j]*lik[i,j]
        
        ubi[(((i-1)*q)+1) : (i*q), i] <- ub
        
        ubbi[(((i-1)*q)+1) : (i*q), (((i-1)*q) + 1) : (i*q)] <- ubb
        
        uybi[(sum(ni[1:i-1])+1) : (sum(ni[1:i])),(((i-1)*q)+1) : (i*q)] <- uyb
        
        uyyi[(sum(ni[1:i-1])+1) : (sum(ni[1:i])),(sum(ni[1:i-1])+1) : (sum(ni[1:i]))] <- uyy
        
        uyi[(sum(ni[1:i-1])+1) : (sum(ni[1:i])),i] <- uy
        
        # Information matrix
        
        # pi
        si.pi <- (Zij[i,j]/as.numeric(pii[j])) - (Zij[i,g]/as.numeric(pii[g]))
        
        # beta
        si.beta <- (1/sigma2[j])*((t(Xi)%*%solve(E)%*%(Zij[i,j]*uy-Zij[i,j]*Ui%*%ub)) - 
                                    (Zij[i,j]*t(Xi)%*%solve(E)%*%Xi%*%beta[,j]))
        #sigmae
        si.sigma <- -(1/2)*(((Zij[i,j]*ni[i])/sigma2[j])-(1/sigma2[j]^2)*((sum(diag(Zij[i,j]*uyy%*%solve(E)))-Zij[i,j]*t(uy)%*%solve(E)%*%mu-Zij[i,j]*t(mu)%*%solve(E)%*%uy-sum(diag(Zij[i,j]*solve(E)%*%((uyb)%*%t(Ui))))-sum(diag(Zij[i,j]*solve(E)%*%((uyb)%*%t(Ui)))) +
                                                                             Zij[i,j]*t(mu)%*%solve(E)%*%Ui%*%ub+Zij[i,j]*t(ub)%*%t(Ui)%*%solve(E)%*%mu+Zij[i,j]*t(mu)%*%solve(E)%*%mu+sum(diag(Zij[i,j]*ubb%*%t(Ui)%*%solve(E)%*%Ui)))))
        
        # phi1 
        Dp <- DevEiAr1(tti,phi1[j],phi2[j],sigma2[j]) 
        Dpr <- Dp$devR_phi1
        dE1_r <- sum(diag(solve(E)%*%Dpr))
        dE2_r <- -(solve(E)%*%Dpr%*%solve(E)) 
        si.phi1 <- - 0.5*Zij[i,j]*dE1_r - (0.5/sigma2[j])*((sum(diag(Zij[i,j]*uyy%*%dE2_r))-Zij[i,j]*t(uy)%*%dE2_r%*%mu-Zij[i,j]*t(mu)%*%dE2_r%*%uy-sum(diag(Zij[i,j]*dE2_r%*%((uyb)%*%t(Ui))))-sum(diag(Zij[i,j]*dE2_r%*%((uyb)%*%t(Ui)))) +
                                                              Zij[i,j]*t(mu)%*%dE2_r%*%Ui%*%ub+Zij[i,j]*t(ub)%*%t(Ui)%*%dE2_r%*%mu+Zij[i,j]*t(mu)%*%dE2_r%*%mu+sum(diag(Zij[i,j]*ubb%*%t(Ui)%*%dE2_r%*%Ui))))      
        # D
        D_der <- DerD(D[[j]])
        deralpha <- rep(0,n_alphas)
        md2<-dim(D[[j]])[1]  
        kont <- 0
        for(i1 in 1:md2){
          for(i2 in 1:(md2+1-i1)){
            kont <- kont+1
            di <- D_der[[i1]][[i2]]
            deralpha[kont] <- (-0.5)*sum(diag(Zij[i,j]*sD[[j]]%*%di-Zij[i,j]*sD[[j]]%*%di%*%sD[[j]]*ubb))     
          }
        }
        
        si <- rbind(si,c(si.beta,si.sigma,si.phi1,deralpha,si.pi))
        
      }# if N
      
      si.g <- cbind(si.g,si)
      
      si <- NULL
      
      ubii[[j]] <- ubi
      beta[,j] <- solve(soma3)%*%(soma4[,j])
      sigma2[j] <- soma2[j]
      D[[j]] <- soma1
      sD[[j]] <- solve(D[[j]]) 
      
      phis <- optimize(f = FCi_phi1,phi22 = phi2[j],
                       beta = beta[,j],sigma2=sigma2[j],tt=tt,ubi=ubi,ubbi=ubbi,uybi=uybi,uyyi=uyyi,uyi=uyi,X=X,U=U,
                       ni=ni,Zij = Zij[,j], interval=c(0.0001,0.9))$minimum
      phi1[j] <- phis
      
      pii[g] <- 1 - (sum(pii) - pii[g])
      
    }# if j
    
    for(v in 1:N) MI = MI + si.g[v,-g*sum(c(p,1,1,1,q*(1+q)/2))]%*%t(si.g[v,-g*sum(c(p,1,1,1,q*(1+ q)/2))])
    
    lk <- sum(log(mixdmvNCens(cc=cc, ni=ni, y=y, tt= tt, X=X, U=U, D=D, beta=beta,phi1= phi1, phi2=phi2, sigma2=sigma2, pii=pii)))
    criteria <- abs((lk/lkante-1))
    lkante <- lk
    
  } #end while
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  
  aic <- -2*lk + 2*d
  bic <- -2*lk + log(N)*d
  edc <- -2*lk + 0.2*sqrt(N)*d
  
  out = list(d=d,sigma2 = sigma2,pii = pii,phi1= phi1,beta = beta,D=D,ubi = ubii,Z = Zij,MI = MI,AIC=aic,BIC=bic,EDC=edc,loglik=lk,time = time.taken)

  return(out)     
  
}

unc_function = function(cc, y, ni, tt, X ,U , g, pii, type_structure, 
                        N, n, p, q, beta, sigma2, D, sD, start.time,
                        error = 0.00001, iter.max = 300){

  cat("Runnig model with", type_structure, "correlation structure ... \n")
  
  theta <- c(beta[,1],sigma2[1],D[[1]][upper.tri(D[[1]], diag = T)],pii[1])
  n_alphas <- length(D[[1]][upper.tri(D[[1]], diag = T)])
  
  d <- g*length(theta) - 1 
  
  start.time <- Sys.time()
  
  criteria <- 1
  count <- 0
  lkante <- 1
  
  while((criteria > error) && (count <= iter.max) ){ # begin while1
  
  si.beta <- si.sigma <- si.pii <- si <- si.g <- NULL
  
  MI <- matrix(0,g*sum(c(p,1,1,q*(1+ q)/2)) - 1 ,g*sum(c(p,1,1,q*(1+ q)/2)) - 1)   
  
  count <- count + 1
  
  print(count)
  
  Zij <- matrix(0, N, g)  
  
  ubii <- list()
  
  ubi <- matrix(0,N*q,N)
  ubbi <- matrix(0,N*q,N*q)
  uybi <- matrix(0,n,N*q)
  uyyi <- matrix(0,n,n)
  uyi <- matrix(0,n,N)
  
  lik <- matrix(0,N,g)
  
  soma5 <- matrix(0,N,g)
  
  for (j in 1:g){
    
    if( sqrt(sum((order(pii)-g:1)*(order(pii)-g:1))) != 0 ){
      betai  <- beta
      sigma2i <- sigma2
      for(l in 1:g){
        beta[,l]    <- betai[,(order(pii, decreasing = TRUE)[l])]
        sigma2[l] <- sigma2i[(order(pii, decreasing = TRUE)[l])]
      }
      pii <- pii[order(pii, decreasing = TRUE)]
    }  
    
    soma1 <- matrix(0,q,q)
    soma2 <- c(rep(0,g)) # It is referred to the sigma2
    soma3 <- matrix(0,p,p)
    soma4 <- matrix(0,p,g) # It is referred to the second sum beta's
    
    
    d1 <- dmvNCens_unc(cc=cc, ni=ni, y=y, tt= tt, X=X, U=U, D=D[[j]], beta=beta[,j], sigma2=sigma2[j])
    if(length(which(d1 == 0)) > 0) d1[which(d1 == 0)] <- .Machine$double.xmin # de donde pertenece
    
    d2 <- mixdmvNCens_unc(cc=cc, ni=ni, y=y, tt= tt, X=X, U=U, D=D, beta=beta, sigma2=sigma2, pii=pii)      
    if(length(which(d2 == 0)) > 0) d2[which(d2 == 0)] <- .Machine$double.xmin # de donde pertenece
    
    Zij[,j] <- (pii[j]*d1)/(d2)
    
    pii[j] <- (1/N)*sum(Zij[,j])
    
    for (i in 1:N){
      
      cci <- cc[(sum(ni[1 :i-1])+1) : (sum(ni[1:i]))]
      yi <- y[(sum(ni[1:i-1])+1) : (sum(ni[1:i]))]
      Xi <- matrix(X[(sum(ni[1:i-1])+1) : (sum(ni[1:i])), ],ncol=p)
      Ui <- matrix(U[(sum(ni[1:i-1])+1) : (sum(ni[1:i])), ],ncol=q)
      tti <- tt[(sum(ni[1:i-1])+1) : (sum(ni[1:i]))]
      mu <- Xi%*%beta[,j]
      Omega <- sigma2[j]*diag(ni[i])
      E <- diag(ni[i])
      Sigma <- (Omega+(Ui)%*%D[[j]]%*%t(Ui))
      Sigma <- (Sigma + t(Sigma))/2
      delta <- solve(sD[[j]]+(t(Ui)%*%solve(E)%*%(Ui*(1/sigma2[j]))))
      delta <- (delta + t(delta))/2
      
      if(sum(cci)==0){
        uy <- matrix(yi,ni[i],1)
        uyy <- yi%*%t(yi)
        ub <- delta%*%(t(Ui)*(1/sigma2[j]))%*%solve(E)%*%(uy-mu)
        ubb <- delta+(delta%*%(t(Ui)*((1/sigma2[j])^2))%*%solve(E)%*%(uyy-uy%*%t(mu) - mu%*%t(uy)+mu%*%t(mu))%*%solve(E)%*%Ui%*%delta)
        uyb <- (uyy-uy%*%t(mu))%*%solve(E)%*%(Ui*(1/sigma2[j]))%*%delta
      }
      
      if(sum(cci)>=1){
        if(sum(cci)==ni[i]){
          muc <- Xi%*%beta[,j]
          Sc <- Sigma
          aux<- Mnormtr(muc,Sc,yi)
          uy<-aux$Ey
          uyy<- aux$Eyy
          ub <- delta%*%(t(Ui)*(1/sigma2[j]))%*%solve(E)%*%(uy-mu)
          ubb <- delta+(delta%*%(t(Ui)*((1/sigma2[j])^2))%*%solve(E)%*%(uyy-uy%*%t(mu) - mu%*%t(uy)+mu%*%t(mu))%*%solve(E)%*%Ui%*%delta)
          uyb <- (uyy-uy%*%t(mu))%*%solve(E)%*%(Ui*(1/sigma2[j]))%*%delta
        }
        else {
          muc <- Xi[cci==1,]%*%beta[,j]+Sigma[cci==1,cci==0]%*%
            solve(Sigma[cci==0,cci==0])%*%(yi[cci==0]-Xi[cci==0,]%*%beta[,j])
          Sc <- Sigma[cci==1,cci==1]-Sigma[cci==1,cci==0]%*%solve(Sigma[cci==0,cci==0])%*%Sigma[cci==0,cci==1]
          Sc <- (Sc + t(Sc))/2
          aux <- Mnormtr(muc,Sc,yi[cci==1])
          uy <- matrix(yi,ni[i],1)
          uy[cci==1] <- aux$Ey
          uyy <- matrix(0,ni[i],ni[i])
          uyy[cci==1,cci==1] <- aux$Vary
          uyy <- uyy+uy%*%t(uy)
          ub <- delta%*%(t(Ui)*(1/sigma2[j]))%*%solve(E)%*%(uy-mu)
          ubb <- delta+(delta%*%(t(Ui)*((1/sigma2[j])^2))%*%solve(E)%*%(uyy-uy%*%t(mu) - mu%*%t(uy)+mu%*%t(mu))%*%solve(E)%*%Ui%*%delta)
          uyb <- (uyy-uy%*%t(mu))%*%solve(E)%*%(Ui*(1/sigma2[j]))%*%delta
        }
      }
      
      soma1 <- soma1 + (Zij[i,j]*ubb)/sum(Zij[,j])          
      
      soma2[j] <- soma2[j] + (sum(Zij[i,j]*diag(uyy%*%solve(E))) - (Zij[i,j]*t(uy)%*%solve(E)%*%mu) - (Zij[i,j]*t(mu)%*%solve(E)%*%uy) - sum(Zij[i,j]*diag(solve(E)%*%((uyb)%*%t(Ui))))-sum(Zij[i,j]*diag(solve(E)%*%((uyb)%*%t(Ui))))
                              +Zij[i,j]*t(mu)%*%solve(E)%*%Ui%*%ub+Zij[i,j]*t(ub)%*%t(Ui)%*%solve(E)%*%mu
                              +Zij[i,j]*t(mu)%*%solve(E)%*%mu+sum(Zij[i,j]*diag(ubb%*%t(Ui)%*%solve(E)%*%Ui)))/(sum(ni[i]*Zij[,j]))
      
      soma3 <- soma3 + Zij[i,j]*(t(Xi)%*%solve(E)%*%Xi)
      
      soma4[,j] <- soma4[,j] + Zij[i,j]*(t(Xi)%*%solve(E)%*%(uy-Ui%*%ub))
      
      soma5[i,j] <- pii[j]*lik[i,j]
      
      ubi[(((i-1)*q)+1) : (i*q), i] <- ub
      
      ubbi[(((i-1)*q)+1) : (i*q), (((i-1)*q) + 1) : (i*q)] <- ubb
      
      uybi[(sum(ni[1:i-1])+1) : (sum(ni[1:i])),(((i-1)*q)+1) : (i*q)] <- uyb
      
      uyyi[(sum(ni[1:i-1])+1) : (sum(ni[1:i])),(sum(ni[1:i-1])+1) : (sum(ni[1:i]))] <- uyy
      
      uyi[(sum(ni[1:i-1])+1) : (sum(ni[1:i])),i] <- uy
      
      # Information matrix
      
      # pi
      si.pi <- (Zij[i,j]/as.numeric(pii[j])) - (Zij[i,g]/as.numeric(pii[g]))
      
      # beta
      si.beta <- (1/sigma2[j])*((t(Xi)%*%solve(E)%*%(Zij[i,j]*uy-Zij[i,j]*Ui%*%ub)) - 
                                  (Zij[i,j]*t(Xi)%*%solve(E)%*%Xi%*%beta[,j]))
      #sigmae
      si.sigma <- -(1/2)*(((Zij[i,j]*ni[i])/sigma2[j])-(1/sigma2[j]^2)*((sum(diag(Zij[i,j]*uyy%*%solve(E)))-Zij[i,j]*t(uy)%*%solve(E)%*%mu-Zij[i,j]*t(mu)%*%solve(E)%*%uy-sum(diag(Zij[i,j]*solve(E)%*%((uyb)%*%t(Ui))))-sum(diag(Zij[i,j]*solve(E)%*%((uyb)%*%t(Ui)))) +
                                                                           Zij[i,j]*t(mu)%*%solve(E)%*%Ui%*%ub+Zij[i,j]*t(ub)%*%t(Ui)%*%solve(E)%*%mu+Zij[i,j]*t(mu)%*%solve(E)%*%mu+sum(diag(Zij[i,j]*ubb%*%t(Ui)%*%solve(E)%*%Ui)))))
      
      # D
      D_der <- DerD(D[[j]])
      deralpha <- rep(0,n_alphas)
      md2<-dim(D[[j]])[1]  
      kont <- 0
      for(i1 in 1:md2){
        for(i2 in 1:(md2+1-i1)){
          kont <- kont+1
          di <- D_der[[i1]][[i2]]
          deralpha[kont] <- (-0.5)*sum(diag(Zij[i,j]*sD[[j]]%*%di-Zij[i,j]*sD[[j]]%*%di%*%sD[[j]]*ubb))     
        }
      }
      
      si <- rbind(si,c(si.beta,si.sigma,deralpha,si.pi))
      
    }
    
    si.g <- cbind(si.g,si)
    
    si <- NULL
    
    ubii[[j]] <- ubi
    beta[,j] <- solve(soma3)%*%(soma4[,j])
    sigma2[j] <- soma2[j]
    D[[j]] <- soma1
    sD[[j]] <- solve(D[[j]]) 
    
    pii[g] <- 1 - (sum(pii) - pii[g])
  }# if j
  
  for(v in 1:N) MI <- MI + si.g[v,-g*sum(c(p,1,1,q*(1+q)/2))]%*%t(si.g[v,-g*sum(c(p,1,1,q*(1+ q)/2))])
  
  lk <- sum(log(mixdmvNCens_unc(cc=cc, ni=ni, y=y, tt= tt, X=X, U=U, D=D, beta=beta, sigma2=sigma2, pii=pii)))
  criteria <- abs((lk/lkante-1))
  lkante <- lk
  
  } #end while
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  
  aic <- -2*lk + 2*d
  bic <- -2*lk + log(N)*d
  edc <- -2*lk + 0.2*sqrt(N)*d
  
  out <- list(d=d, sigma2 = sigma2, pii = pii, beta = beta, D=D, ubi = ubii,Z = Zij,MI = MI,AIC=aic,BIC=bic,EDC=edc,loglik=lk,time = time.taken)
  
  return(out)     
  
}

## List with correlation structure functions ##
funclist = list("dec" = dec_function, "Phi1" = Phi1_function, "unc" = unc_function)




