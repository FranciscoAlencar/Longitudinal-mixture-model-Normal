#--------------------------------------------#
# EM algorithm for the FINITE MIXTURE LINEAR #
# MIXED MODEL FOR CENSORED DATA              #
# Creatred by Alencar, Matos and Lachos      #
#--------------------------------------------#

#@----------------------#@
#@  AR - Simulation     #@
#@----------------------#@

tN <- 200 # Number of subjects
tni <- c(rep(8,tN)) # Number of observations for each subject yi with the same number.
tn <- sum(tni)

tpii <- c(0.7,0.3) # pii

tg <- length(tpii) # Number of mixture elements

tpii <- tpii[order(tpii, decreasing = TRUE)]

tphi1 <- c(0.6,0.4)

tphi2 <- c(1,1)

tsigma2 <- c(1,2)

tbeta1 <- matrix(c(-3,-1,-4,-2),ncol = 1) 

tp <- length(tbeta1) # Number of betas

tbeta2 <- matrix(c(1,2,3,4),ncol = 1) 

beta = matrix(0, ncol = 2, nrow = tp)

beta[,1] = tbeta1
beta[,2] = tbeta2

ttt <- matrix(rep(c(0,1,3,6,9,12,18,24),tN),ncol =1)

tD1 <- matrix(c(0.5,0.1,0.1,0.5), ncol = 2)

tD2 <- matrix(c(0.5,0.1,0.1,0.5), ncol = 2)

D = list(tD1, tD2)

tX <- matrix(c(rep(1,sum(tni)),runif(tN*tni[1],1,5),runif(tN*tni[1],0,1),
               runif(tN*tni[1],-1,1)),ncol = tp)

tU <- tX[,1:2]

tq <- dim(tU)[2] # Number of b's

#seed <- 6584

pCens <- c(0.05, 0.05)

N_elements = tN; n_obs = tni;  X_matrix = tX;  U_matrix = tU; beta_matrix = beta; 
pi_vec = tpii; phi1_vec = tphi1; phi2_vec = tphi2; sigma2_vec = tsigma2;
time_vec = ttt; D_matrix_list = D; pCens_vec = pCens
#seed_rep = seed; each_group = FALSE 
