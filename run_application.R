rm(list=ls(all=TRUE))

setwd("~/Documentos/Hildemar_notebook/Documents/Doutorado/tese/R/FM_LMM")
source("data_application.R")
source("fit_Nfmlmm.R")
source("load_packages.R")

AA_VH <- Cens.fmlmEM(cc = tcc, y = ty, ni= tni, tt= ttt, type_structure = "unc",
                     X = tX ,U = tU , g = 1, error = 0.00001, iter.max = 350)

AA_VH$AIC
AA_VH$BIC
AA_VH$EDC
