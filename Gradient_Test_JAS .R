#Program Example 3 -JAS
# Idemauro Antonio Rodrigues de Lara - Cesar Augusto Taconeli
#With 4 time occasions

rm(list=ls(all=TRUE))

#######################################################################################################
# the data-set

setwd("~/Documents/Pesquisa/Orientandos/Concluídos/Doutorado/Laura/Artigo2_JAS/Terceira_submissao/exemplo3")
#setwd("C:/Users/Idemauro/Documents/Pesquisa/Exemplo3_JAS")
dados <- read.csv("dados_adaptados.csv", head=TRUE, sep=";", dec=",")
dados=na.omit(dados)
head(dados)
dim(dados)

dados <- dados[,c('instar','planta','resp1','resp2','resp3','resp4')]
dados$id <- seq(1,nrow(dados))

# Conversion for factors

dados$instar <- as.factor(dados$instar)
dados$planta <- as.factor(dados$planta)
dados$resp1 <- as.factor(dados$resp1)
dados$resp2 <- as.factor(dados$resp2)
dados$resp3 <- as.factor(dados$resp3)
dados$resp4 <- as.factor(dados$resp4)

# Packages
require(MASS)
require(car)
require(VGAM)
require(Matrix)
require(nnet)
require(ordinal)

#######################################################################################################
# Part 2 - Functions to change structure of data

widetolong<-function(basewide,varswide){
  r1<-reshape(basewide,varying=list(varswide),direction='long',v.names='resp')
  r1$time<-factor(r1$time)
  rownames(r1)<-1:nrow(r1)
  r1<-r1[,-which(names(r1)=='id')]
  r2<-r1[-which(r1$time==levels(r1$time)[1]),]
  r2<-cbind(r2,r1[-which(r1$time==levels(r1$time)[length(levels(r1$time))]),'resp'])
  names(r2)[ncol(r2)]<-'respp';head(r2)
  return(r2)
}

longtowide<-function(baselonga,resplonga,vartime,respprevia,idsuj,resp1){
  indicat1<-names(which(table(baselonga$time)!=0))[1]
  r3<-cbind(reshape(baselonga[,-which(names(baselonga)==respprevia)],idvar = idsuj,v.names = resplonga,
                    timevar = vartime, direction = "wide",sep=""),
            baselonga[,respprevia][which(baselonga[,vartime]==indicat1)])
  names(r3)[ncol(r3)]<-resp1
  return(r3)
}

#######################################################################################################
dadosl <- widetolong(dados,c('resp1','resp2','resp3','resp4')) # data in long way.

dadosl$resp <- factor(dadosl$resp)
dadosl$respp <- factor(dadosl$respp)
dadosl$resp=relevel(dadosl$resp,ref="nothing")
dadosl$respp=relevel(dadosl$respp,ref="nothing")
str(dadosl)
dados$resp1=relevel(dados$resp1,ref="nothing")
dados$resp2=relevel(dados$resp2,ref="nothing")
dados$resp3=relevel(dados$resp3,ref="nothing")
dados$resp4=relevel(dados$resp4,ref="nothing")

################################################################################################
# Proposed Gradient test - using clm
################################################################################################
#Models
modEg <- clm(resp ~ 1, nominal= ~ planta+respp, data=dadosl)

mod12g<-clm(resp2~ 1, nominal= ~planta+resp1,data=dados)

mod23g<-clm(resp3~ 1, nominal= ~planta+resp2,data=dados)

mod34g<-clm(resp4~ 1, nominal= ~planta+resp3,data=dados)


#Coeficients

beta0g <- coef(modEg)
beta12g <- coef(mod12g)
beta23g <- coef(mod23g)
beta34g <- coef(mod34g)

beta <- c(beta12g,beta23g,beta34g)
beta0vec <- rep(beta0g,3)

length(beta0g)

############################ Gradient statistic ################################


getFittedC <- ordinal:::getFittedC

clm.nll <- function(rho, par) {
  if(!missing(par)) rho$par <- par
  with(rho, {
    if(k > 0)
      sigma <- Soff * exp(drop(S %*% par[n.psi + 1:k]))
    ### NOTE: we have to divide by sigma even if k=0 since there may be an
    ### offset but no predictors in the scale model:
    eta1 <- (drop(B1 %*% par[1:n.psi]) + o1)/sigma
    eta2 <- (drop(B2 %*% par[1:n.psi]) + o2)/sigma
  })
  ### NOTE: getFitted is not found from within rho, so we have to
  ### evalueate it outside of rho
  rho$fitted <- getFittedC(rho$eta1, rho$eta2, rho$link, rho$par[length(rho$par)])
  if(all(is.finite(rho$fitted)) && all(rho$fitted > 0))
    ### NOTE: Need test here because some fitted <= 0 if thresholds are
    ### not ordered increasingly.
    -sum(rho$wts * log(rho$fitted))
  else Inf
}

clm.grad <- function(rho) {
  ### requires that clm.nll has been called prior to
  ### clm.grad.
  with(rho, {
    p1 <- if(!nlambda) dfun(eta1) else dfun(eta1, lambda)
    p2 <- if(!nlambda) dfun(eta2) else dfun(eta2, lambda)
    wtpr <- wts/fitted
    C2 <- B1*p1/sigma - B2*p2/sigma
    if(k <= 0) return(-crossprod(C2, wtpr))
    C3 <- -(eta1 * p1 - eta2 * p2) * S
    return(-crossprod(cbind(C2, C3), wtpr))
    ### NOTE: C2 and C3 are used by clm.hess
  })
}

clm.grad_direct <- function(rho, par) {
  ### does not require that clm.nll has been called prior to
  ### clm.grad.
  clm.nll(rho, par)
  clm.grad(rho)
}

################################################################################
modEe<- update(modEg, doFit=FALSE)
estimatE <- coef(modEg) 


mod12a <- update(mod12g, doFit=FALSE)
mod23a <- update(mod23g, doFit=FALSE)
mod34a <- update(mod34g, doFit=FALSE)


### Gradient function using liklihood estimates by clm

clm.grad_direct(modEe, estimatE)
modEg$gradient

vetescore <- c(clm.grad_direct(mod12a, estimatE),
               clm.grad_direct(mod23a, estimatE),
               clm.grad_direct(mod34a, estimatE))
               

grad <- abs(vetescore%*%(beta-beta0vec));round(grad,4)
pchisq(grad, 40-8, lower.tail = FALSE)


