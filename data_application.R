#@----------------------#@
#@   Application-AR_VH  #@
#@----------------------#@

data1 <- read.csv("dataA5055.csv")

attach(data1)

data1 <- subset(data1, !is.na(cd4))
data1['cd4']

length(cd8[!is.na(cd8)])

data1 <- subset(data1, !is.na(cd8))
data1['cd8']

subjects <- unique(data1$patid)
length(subjects)
table(data1['patid'])

cluster <- c(match(data1$patid,subjects))
tN <- length(subjects) # quantidade de indivíduos
n <-length(cluster) # total de observações
y1 <- c(data1$logrna) # Nossa variável resposta
y2.1 <- c(data1$cd4)
y2.2 <- c(data1$cd8)

y2 <- y2.1/y2.2
y2<-log(y2)

tni <- c()

for (j in 1:tN) {
  tni[j] <- sum(cluster==j)
}

new_time_aux <- c(0, 7, 14, 28, 56, 84, 112, 140, 168)

new_time <- NULL

for(i in 1:tN){
  new_time <- c(new_time,new_time_aux[1:tni[i]])  
}

# Dados antes da censura
datas_sc <- cbind(data1$patid,as.matrix(cluster,tN,1),data1$arm,data1$day,
                  as.matrix(y1,tN,1),as.matrix(y2,tN,1),data1$cd4,data1$cd8)
nam_row=as.character((1:n))
nam_col=c("patid","cluster","arm","day","logrna","logcd48","cd4","cd8")

datas_sc <- matrix(datas_sc,nrow=n,ncol=8,dimnames=list(" "=nam_row," "=nam_col))
dados_sc <- as.data.frame(datas_sc)
attach(dados_sc)

dados2_sc <- dados_sc

dados2_sc$arm[dados2_sc$arm=="1"] <- "treatment 1"
dados2_sc$arm[dados2_sc$arm=="2"]   <- "treatment 2"
head(dados2_sc, 3)

cc <- (data1$rna<50)+0
y1[y1<=log10(50)] <- log10(50) # Censura

datas <- cbind(data1$patid,as.matrix(cluster,tN,1),data1$arm,data1$day,
               as.matrix(y1,tN,1),as.matrix(y2,tN,1),as.matrix(cc,tN,1),data1$cd4)
nam_row=as.character((1:n))
nam_col=c("patid","cluster","arm","day","logrna","logcd48","cens","cd4")

datas <- matrix(datas,nrow=n,ncol=8,dimnames=list(" "=nam_row," "=nam_col))
dados <- as.data.frame(datas)
attach(dados)

dados2 <- dados
dados2$arm[dados2$arm=="1"] <- "treatment 1"
dados2$arm[dados2$arm=="2"]   <- "treatment 2"
head(dados2, 3)

# Mean per treatment 
dados_treat1 = subset(dados2, dados2$arm== "treatment 1")
aux1_treat1  = list()
n_treat1 = length(unique(dados_treat1$patid))
nt1 = NULL
for (j in unique(dados_treat1$cluster)) {
  nt1 <- c(nt1,sum(dados_treat1$cluster==j))
}
Cont1 = 0
for(i in unique(dados_treat1$cluster)){
  Cont1 = Cont1 + 1
  aux1_treat1[[i]] = t(dados_treat1$logrna[dados_treat1$cluster == i])
}
aux1_treat1 

library(plyr)
M_treat1 <- t(rbind.fill.matrix(aux1_treat1))
MM_treat1 <- c()
for(i in 1:9){
  MM_treat1[i] <-  mean(M_treat1[i,], na.rm = TRUE) 
}


dados_treat2 = subset(dados2, dados2$arm== "treatment 2")
aux1_treat2  = list()
n_treat2 = length(unique(dados_treat2$patid))
nt2 = NULL
for (j in unique(dados_treat2$cluster)) {
  nt2 <- c(nt2,sum(dados_treat2$cluster==j))
}
Cont2 = 0
for(i in unique(dados_treat2$cluster)){
  Cont2 = Cont2 + 1
  aux1_treat2[[Cont2]] = t(dados_treat2$logrna[dados_treat2$cluster == i])
}
aux1_treat2

M_treat2 <- t(rbind.fill.matrix(aux1_treat2))
MM_treat2 <- c()
for(i in 1:9){
  MM_treat2[i] <-  mean(M_treat2[i,], na.rm = TRUE) 
}

# Censuring per treatments
to_treat1 <- sum(dados$arm == 1)
Cens_treat1 <- round(sum(subset(dados$arm == 1, dados$cens==1))/to_treat1,3)*100

to_treat2 <- sum(dados$arm == 2)
Cens_treat2 <- round(sum(subset(dados$arm == 2, dados$cens==1))/to_treat2,3)*100

dados_treat1 = subset(dados2, dados2$arm== "treatment 1")
dados_treat1$arm[dados_treat1$arm == "treatment 1"] <- "Treatment 1"

TIME = c(0,   7,  28,  56,  89, 105, 141, 168, 197)

MM_treat1_data = data.frame(y = MM_treat1, x = TIME)

dados_treat2 = subset(dados2, dados2$arm== "treatment 2")
dados_treat2$arm[dados_treat2$arm == "treatment 2"] <- "Treatment 2"

MM_treat2_data = data.frame(y = MM_treat2, x = TIME)

aux_data = rep("Treatments",dim(dados)[1])
dados2_new <- cbind(dados2,aux_data)

# AJUSTE 

tcc <- c(dados2$cens)
ty <- dados2$logrna
ttt <- dados2$day/7

new_cd4 = (data1$cd4 - mean(data1$cd4))/sd(data1$cd4)

# model 1
tX <- matrix(c(rep(1,length(ty)),ttt,ttt^0.5),nrow = length(ty),ncol = 3)
dim(tX)
#tU <- matrix(c(rep(1,length(ty)),ttt),nrow = length(ty),ncol = 2)

tU <- matrix(c(rep(1,length(ty))),nrow = length(ty),ncol = 1)
