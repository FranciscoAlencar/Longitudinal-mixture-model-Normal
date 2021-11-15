rm(list=ls(all=TRUE))

setwd("~/Documentos/Hildemar_notebook/Documents/Doutorado/tese/R/FM_LMM")
source("simulated_data.R")
source("simulation_function.R")
source("fit_Nfmlmm.R")
source("load_packages.R")

seed = 6584

seed_rep = set.seed(seed)

Data_simu = gen_norm_long_cens(N_elements, n_obs, X_matrix, U_matrix, beta_matrix, 
                   pi_vec, phi1_vec, phi2_vec, sigma2_vec, time_vec, 
                   D_matrix_list, pCens_vec, seed_rep, each_group = TRUE)

ty = Data_simu$y
tcc = Data_simu$censured

hist(ty, freq = TRUE, breaks = 30)
  
AA  = Cens.fmlmEM(cc = tcc, y = ty, ni= n_obs, tt= time_vec, type_structure = "ar",
X = X_matrix ,U = U_matrix , g = 2, error = 0.0001, iter.max = 300)
  