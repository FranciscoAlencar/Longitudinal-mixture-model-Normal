##-------------------------------------------------------##
##      Generate censored longitudinal data function     ##
##-------------------------------------------------------##

gen_norm_long_cens <- function(N_elements, n_obs, X_matrix, U_matrix, beta_matrix, 
                   pi_vec, phi1_vec, phi2_vec, sigma2_vec, time_vec, 
                   D_matrix_list, pCens_vec, seed_rep, each_group){
  
  tg = length(pi_vec)  
  tp = dim(beta_matrix)[1] # Number of betas
  tq = dim(U_matrix)[2] # Number of b's
  tn = sum(n_obs)
  
  Y = list()
  CC = list()
  
  for(j in 1:tg){
    Y[[j]] = c(rep(NaN,sum(n_obs)))
    CC[[j]] = c(rep(NaN,sum(n_obs)))
  }
  
  Y[["y"]] = c(rep(rep(0,sum(n_obs))))
  CC[["cc"]] = c(rep(rep(0,sum(n_obs))))

  set.seed(seed_rep)
  P = rmultinom(N_elements, size = 1, prob = pi_vec)
  
  if(each_group){
  
    for(i in 1:N_elements){
      
      group = which.max(P[,i])
      
      set.seed(seed_rep)
      bi =  matrix(rmvnorm(1, mean = rep(0,tq), D_matrix_list[[group]]), ncol = 1)
      
      tOmega = DecM(time_vec[(sum(n_obs[1:i-1])+1) : (sum(n_obs[1:i]))], phi1_vec[group], phi2_vec[group], sigma2_vec[group])
       
      set.seed(seed_rep)
      e = matrix(rmvnorm(1,mean = rep(0, n_obs[i]), tOmega), ncol =1)
      
      Y[[group]][(sum(n_obs[1:i-1])+1) : (sum(n_obs[1:i]))] = matrix(X_matrix[(sum(n_obs[1:i-1])+1) : (sum(n_obs[1:i])), ],ncol = tp)%*%beta_matrix[,group] + 
        matrix(U_matrix[(sum(n_obs[1:i-1])+1) : (sum(n_obs[1:i])), ], ncol=tq)%*%bi + e
        
    }
    
    for(j in 1:tg){
      aa <- sort(Y[[j]],decreasing = FALSE)
      cutof <- aa[ceiling(pCens_vec[j]*length(Y[[j]]))]
      CC[[j]] <- matrix(1,length(Y[[j]]),1)*(Y[[j]] <= cutof)
      Y[[j]][CC[[j]]==1] = cutof
      CC[[1]][which(is.na(CC[[1]]))] <- CC[[j]][which(is.na(CC[[1]]))]
    }
    
    CC[["cc"]] = CC[[1]]
    
    for(u in 1:N_elements){
      
      group = which.max(P[,u])
      
      Y[["y"]][(sum(n_obs[1:u-1])+1) : (sum(n_obs[1:u]))] = Y[[group]][(sum(n_obs[1:u-1])+1) : (sum(n_obs[1:u]))]
      
    }

  }else{
    
    for(i in 1:N_elements){
      
      group = which.max(P[,i])
      
      bi =  matrix(rmvnorm(1, mean = rep(0,tq), D_matrix_list[[group]]), ncol = 1)
      
      tOmega = DecM(time_vec[(sum(n_obs[1:i-1])+1) : (sum(n_obs[1:i]))], phi1_vec[group], phi2_vec[group], sigma2_vec[group])
      
      e = matrix(rmvnorm(1,mean = rep(0, n_obs[i]), tOmega), ncol =1)
      
      Y[["y"]][(sum(n_obs[1:i-1])+1) : (sum(n_obs[1:i]))] = matrix(X_matrix[(sum(n_obs[1:i-1])+1) : (sum(n_obs[1:i])), ],ncol = tp)%*%beta_matrix[,group] + 
        matrix(U_matrix[(sum(n_obs[1:i-1])+1) : (sum(n_obs[1:i])), ], ncol=tq)%*%bi + e
      
    }
    
    aa = sort(Y[["y"]],decreasing = FALSE)
    cutof = aa[ceiling(pCens_vec[1]*N_elements*n_obs[1])]
    CC[["cc"]] = matrix(1,length(Y[["y"]]),1)*(Y[["y"]] <= cutof)
    
    Y[["y"]][CC[["cc"]]==1] = cutof
    
  }
 
  out = list(y = Y[["y"]], censured = CC[["cc"]])
  
  return(out)
      
}

  