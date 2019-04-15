#-----------------------------------------------------------
#
# NRPois-fun.R
#
# Function for Newton-Raphson method for Poisson data
#
# Programmed by Kasun Yang & Sang Kyu Lee at 2019-01-03
#
#-----------------------------------------------------------
nr.pois <- function(xij, wij, theta0, yi, alpha0, select.dist, eps, itmax, trunc.num) {
  n <- nrow(xij); pdim <- ncol(xij); qdim <- ncol(wij)
  n1 <- n-sum(is.nan(yi))
  if(n<=n1) return(list(theta<-NA, Jbar<-NA, likval<-NA))
  pardim <- pdim+qdim
  
  if(select.dist=="gumbel") {
    fga <- function(ga02) { 
      temp1 <- sum(pexp( exp(ga02%*%t(wij[1:n1,])) , log.p=TRUE))
      temp2 <- -sum(exp(  ga02%*%t(wij[(n1+1):n,]) ))
      return(loglik <- -(temp1+temp2))
    }
  } else if(select.dist=="normal") {
    fga <- function(ga02) {
      temp1 <- sum(pnorm( ga02%*%t(wij[1:n1,]) ,log.p=TRUE))
      temp2 <- sum(pnorm( ga02%*%t(wij[(n1+1):n,]) , lower.tail=FALSE, log.p=TRUE))
      return(loglik <- -(temp1+temp2))
    }
  } else if (select.dist=="logistic") {
    fga <- function(ga02) {
      temp1 <- sum(plogis( ga02%*%t(wij[1:n1,]) ,log.p=TRUE))
      temp2 <- sum(plogis( ga02%*%t(wij[(n1+1):n,]) , lower.tail=FALSE, log.p=TRUE))
      return(loglik <- -(temp1+temp2)) # addition done
    }
  }
  theta <- theta0
  
  if (is.null(theta0)) # To give initial values
  {
    observed <- !is.na(yi)
    b0b1glm <- glm(yi[1:n1] ~ xij[1:n1,2:pdim], family=poisson(log))
    ga_esti <- glm(as.numeric(observed) ~ wij[,2:qdim], family=binomial(probit))
#      if(select.dist=="gumbel") {
#        ga_esti <- glm(as.numeric(observed) ~ wij[,2:qdim], family=binomial(logit))
#      } else {
#        ga_esti <- glm(as.numeric(observed) ~ wij[,2:qdim], family=binomial(probit))
#      }
#    ga_esti <- optim(c(rep(1,qdim)),fga)$par
#    theta <- as.numeric(c(b0b1glm$coef,ga_esti)) # initial values)
    theta <- as.numeric(c(b0b1glm$coef,ga_esti$coef)) # initial values)
  }

  alpha <- alpha0
  diff <- 1
  it <- 0; initmax <- 10000 # To avoid an infinite loop
  likval <- like_poi(xij, wij, theta, yi, alpha, select.dist, trunc.num)
  if(is.na(likval)) {
    return(list(theta<-NA, Jbar<-NA, likval<-NA)) # Key Point where NA occurs #
  }
  
  while ((diff > eps) && (it < itmax))
  { 
    it <- it+1
    theta.old <- theta; likval.old <- likval
    
    # score and observed information matrix ------------------------------------------------------
    score_Jbar <- eval_score_obsinf_poi(xij, wij, theta, yi, alpha, select.dist, trunc.num) 
    
    score <- score_Jbar[[1]]
    Jbar <- score_Jbar[[2]]
    
    tCh <- tryCatch(is.positive.definite(Jbar),error=function(e){999})
    if(tCh != 999) {
      Jbar <- posdefify(Jbar) # To stabilize additionally (make it a pd matrix)
    } else {
      return(list(theta<-NA, Jbar<-NA, likval<-NA)) # Key Point where NA occurs #
    }

    likval <- likval.old - 1; delta <- 1
    init <- 0
    while((likval < likval.old) && (init < initmax)) {
      
      theta <- theta.old + delta*solve(Jbar,score)
      init <- init+1
      
      #       ifelse(alpha < alpha0,likval <- -Inf,likval <- like(mi, xij, wij, theta, yi, alpha, negasso))  
      
      likval <- like_poi(xij, wij, theta, yi, alpha, select.dist, trunc.num)

      if(is.na(likval)) likval <- -Inf
      
      delta <- delta/2
    } # end of while((likval < likval.old) && (init < initmax))
    diff <- sum(abs(theta-theta.old)) 
#cat("within while",it,likval,"\n")
  } # end of while ((diff > eps) && (it < itmax))
  list(theta<-theta, Jbar<-Jbar, likval<-likval) # estimated parameters and Fisher observed information matrix
} # end of nr.pois, done

like_poi <- function(xij, wij, theta, yi, alpha, select.dist, trunc.num) {
  # log-likelihood function ---------------------------------------------
  pdim <- ncol(xij); qdim <- ncol(wij)
  pardim <- pdim+qdim
  n <- nrow(xij)
  n1 <- n-sum(is.nan(yi)) # count # of NaN's
  
  xibeta <- colSums(diag(theta[1:pdim])%*%t(xij))
  wigamma <- colSums(diag(theta[(pdim+1):(pdim+qdim)])%*%t(wij))
  mui <- exp(xibeta)
  
  alpha <- alpha
  taui <- wigamma
  etai <- alpha*exp(-xibeta)
#  etai <- alpha/mui
  k <- 0:trunc.num
  
  if(select.dist=="gumbel") {
    hyi <- exp(taui+etai*yi) 
  } else {
    hyi <- taui+etai*yi
  }
  
  tempLik1 <- c()
  
  if(select.dist=="gumbel"){
    tempLik0 <- dpois(yi[1:n1], mui[1:n1], log=TRUE) + pexp(hyi[1:n1], log.p=TRUE)
    for(i in (n1+1):n) {
      tempLik1[i-n1] <- log( sum( exp(dpois(k,mui[i],log=TRUE))*exp(-exp(taui[i]+etai[i]*k)) ) + .Machine$double.eps)
    }
  } else if (select.dist=="normal") {
    tempLik0 <- dpois(yi[1:n1], mui[1:n1], log=TRUE) + pnorm(hyi[1:n1], log.p=TRUE)
    for(i in (n1+1):n) {
      tempLik1[i-n1] <- log(1-sum( exp(dpois(k,mui[i], log=TRUE))*exp(pnorm(taui[i]+etai[i]*k, log.p=TRUE)) )+.Machine$double.eps)
    }
  } else if (select.dist=="logistic"){
    tempLik0 <- dpois(yi[1:n1], mui[1:n1], log=TRUE) + plogis(hyi[1:n1], log.p=TRUE)
    for(i in (n1+1):n) {
      tempLik1[i-n1] <- log(1-sum( exp(dpois(k,mui[i], log=TRUE))*exp(plogis(taui[i]+etai[i]*k, log.p=TRUE)) )+.Machine$double.eps)
    }
  }
  tempLik <-  c(tempLik0, tempLik1)  
  
  return(sum(tempLik))
} # done

eval_score_obsinf_poi <- function(xij, wij, theta, yi, alpha, select.dist, trunc.num) {
  # score function ---------------------------------------------------
  n <- nrow(xij); pdim <- ncol(xij); qdim <- ncol(wij)
  pardim <- pdim+qdim 
  n1 <- n-sum(is.nan(yi)) # count # of NaN's
  
  score <- c()
  pii <- c()
  ppipmui <- c()
  ppiptaui <- c()

  xibeta <- colSums(diag(theta[1:pdim])%*%t(xij))
  wigamma <- colSums(diag(theta[(pdim+1):(pdim+qdim)])%*%t(wij))
  mui <- exp(xibeta)
  
  alpha <- alpha
  taui <- wigamma
  etai <- alpha*exp(-xibeta)
#  etai <- alpha/mui
  k <- 0:trunc.num
  
  if(select.dist=="gumbel") {
    hyi <- exp(taui+etai*yi) 
    g0G0 <- exp( dexp(hyi[1:n1],log=TRUE) - pexp(hyi[1:n1], log.p=TRUE) )
    g0pG0 <- -g0G0
    phyipmui <- -etai[1:n1]*yi[1:n1]*exp(-xibeta[1:n1])*hyi[1:n1]
    phyiptaui <- hyi[1:n1]
  } else if(select.dist=="normal") {
    hyi <- taui+etai*yi
    # First calculate log values and then exponentiate    
    g0G0 <- exp( dnorm(hyi[1:n1],log=TRUE) - pnorm(hyi[1:n1],log.p=TRUE) )
    g0pG0 <- -hyi[1:n1]*g0G0
    phyipmui <- -etai[1:n1]*yi[1:n1]*exp(-xibeta[1:n1])
    phyiptaui <- rep(1,n1)
  } else if(select.dist=="logistic") {
    hyi <- taui+etai*yi
    g0G0 <- exp( dlogis(hyi[1:n1],log=TRUE)-plogis(hyi[1:n1],log.p=TRUE) ) 
    g0pG0 <- (1/(1+exp(hyi[1:n1])))*((1-exp(hyi[1:n1]))/(1+exp(hyi[1:n1])))
    phyipmui <- -etai[1:n1]*yi[1:n1]/mui[1:n1]
    phyiptaui <- rep(1,n1)
  }
  
  if (select.dist=="gumbel") {
    for(i in (n1+1):n) {
      tek <- taui[i]+etai[i]*k
      edpk <- exp(dpois(k,mui[i],log=TRUE))
      pii[i-n1] <- 1- sum( edpk*exp(-exp(tek)) ) #- .Machine$double.eps # Note
      if (pii[i-n1]==1) pii[i-n1] <- 1 - .Machine$double.eps # to avoid 1
      ppipmui[i-n1] <- sum( edpk*exp(-exp(tek))*(1-k/mui[i]-etai[i]*k/mui[i]*exp(tek)) )
      ppiptaui[i-n1] <- sum( edpk*exp(tek)*exp(-exp(tek)) )
    }
  } else if(select.dist=="normal") {
    for (i in (n1+1):n) {
      tek <- taui[i]+etai[i]*k
      edpk <- exp(dpois(k,mui[i],log=TRUE))
      eptek <- exp(pnorm(tek,log.p=TRUE))
      edtek <- exp(dnorm(tek,log=TRUE))
      pii[i-n1] <- sum( edpk*eptek ) #- .Machine$double.eps  # Note
      if (pii[i-n1]==1) pii[i-n1] <- 1 - .Machine$double.eps # to avoid 1 # Note
      ppipmui[i-n1] <- sum( edpk*( (k/mui[i]-1)*eptek-etai[i]*k/mui[i]*edtek ) )
      ppiptaui[i-n1] <- sum( edpk*edtek )
    }
  } else if (select.dist=="logistic"){
    for(i in (n1+1):n){
      pii[i-n1] <- sum(exp(dpois(k,mui[i],log=TRUE) + plogis(taui[i]+etai[i]*k,log.p=TRUE ) ) )
      
      ppipmui[i-n1] <- sum( exp( dpois(k, mui[i],log = TRUE) + plogis(taui[i]+etai[i]*k, log.p = TRUE) ) * (k/mui[i] -1 -
                                                                                                           plogis(taui[i]+etai[i]*k)*exp(-(taui[i]+etai[i]*k) )*(etai[i]*k/mui[i]) )   )
      ppiptaui[i-n1] <- sum( exp(dpois(k, mui[i],log=TRUE)) *  (plogis(taui[i]+etai[i]*k))^2*
                            exp(-(taui[i]+etai[i]*k)))
    }
  }
  
  for (j in 1:pdim) {
    score[j] <- t( ( (yi[1:n1]-mui[1:n1])/mui[1:n1] + g0G0*phyipmui)*mui[1:n1] )%*%xij[1:n1,j] +
      t( (ppipmui/(pii-1))* mui[(n1+1):n] )%*%xij[(n1+1):n,j]
  }
  
  for (j in 1:qdim) {
    
    if ( (sum(g0G0==0))==n1 ) {
      score[pdim+j] <- t(ppiptaui/(pii-1))%*%wij[(n1+1):n,j] 
    } else {
      score[pdim+j] <- t(g0G0*phyiptaui)%*%wij[1:n1,j] + t(ppiptaui/(pii-1))%*%wij[(n1+1):n,j] 
    }  
  }
  
  # observed information matrix  ------------------------------------------------------
  Jbar <- matrix(0,pardim,pardim)
  p2pipmui2 <- c()
  p2piptaumui <- c()
  p2piptaui2 <- c()
  
  
  if (select.dist=="gumbel") {
    p2hyipmui2 <- (etai[1:n1]*yi[1:n1]*exp(-2*xibeta[1:n1]))*hyi[1:n1]*(2+etai[1:n1]*yi[1:n1])
    p2hyiptaumui <- -etai[1:n1]*yi[1:n1]*hyi[1:n1]*exp(-xibeta[1:n1])
    p2hyiptaui2 <- hyi[1:n1] 
  } else if (select.dist=="normal"){
    p2hyipmui2 <- 2*etai[1:n1]*yi[1:n1]*exp(-2*xibeta[1:n1])
    p2hyiptaumui <- rep(0,n1)
    p2hyiptaui2 <- rep(0,n1)
  } else if (select.dist=="logistic") {
    p2hyipmui2 <- 2*yi[1:n1]*etai[1:n1]/mui[1:n1]^2
    p2hyiptaumui <- rep(0,n1)
    p2hyiptaui2 <- rep(0,n1)
  }
  
  if (select.dist=="gumbel") {
    for(i in (n1+1):n) {
      edpk2 <- exp(dpois(k,mui[i],log=TRUE))
      etek2 <- exp(taui[i]+etai[i]*k)
      p2pipmui2[i-n1] <- sum( edpk2*exp(-etek2)*((1-k/mui[i]-etai[i]*k/mui[i]*etek2)^2 +k/mui[i]^2*(1+etai[i]*etek2*(2+etai[i]*k))))
      p2piptaumui[i-n1] <- sum( edpk2*etek2*exp(-etek2)*(k/mui[i]-1-(etai[i]*k/mui[i])*(1-etek2)))      
      p2piptaui2[i-n1] <- sum( edpk2*etek2*exp(-etek2)*(1-etek2))
    }
  } else if (select.dist=="normal") {
    for(i in (n1+1):n) {
      edpk2 <- exp(dpois(k,mui[i],log=TRUE))
      tek2 <- taui[i]+etai[i]*k
      edtek2 <- exp(dnorm(tek2,log=TRUE))
      eptek2 <- exp(pnorm(tek2,log.p=TRUE))
      p2pipmui2[i-n1] <- sum( edpk2*( eptek2*((k/mui[i]-1)^2-k/mui[i]^2) - etai[i]*k/mui[i]*edtek2*( 2*k/mui[i]-2-2/mui[i]+etai[i]*k/mui[i]*tek2 ) ))
      p2piptaumui[i-n1] <- sum( edpk2*edtek2*( k/mui[i]-1+etai[i]*k/mui[i]*tek2 ))
      p2piptaui2[i-n1] <- - sum( edpk2*tek2*edtek2)
    }
  } else if (select.dist=="logistic") {
    for (i in (n1+1):n) {
      star <- plogis(taui[i]+etai[i]*k)*exp(-(taui[i]+etai[i]*k))*(etai[i]*k/mui[i])
      p2pipmui2[i-n1] <- sum(  exp(dpois(k,mui[i],log=TRUE) + plogis(taui[i]+etai[i]*k,log.p=TRUE ) )*( (k/mui[i]-1)*(k/mui[i]-1-star)-star*(k/mui[i]-1-star)+(-k/mui[i]^2+star*(star-etai[i]*k/mui[i]+2/mui[i])) )  )
      p2piptaumui[i-n1] <- sum( dpois(k,mui[i])*((plogis(taui[i]+etai[i]*k))^2*
                                                exp(-(taui[i]+etai[i]*k))*(k/mui[i] -1 -star)
                                              -(plogis(taui[i]+etai[i]*k))^2*exp(-(taui[i]+etai[i]*k))*(star-etai[i]*k/mui[i])) )
      p2piptaui2[i-n1] <- sum( dpois(k,mui[i]) * ((2*(plogis(taui[i]+etai[i]*k))^3)*
                                                 (exp(-2*(taui[i]+etai[i]*k)))-(plogis(taui[i]+etai[i]*k))^2*exp(-(taui[i]+etai[i]*k)) ))
    }
  }

  p2likb0b0Temp0 <- -1+(g0pG0-g0G0^2)*phyipmui^2*mui[1:n1]+ g0G0*( p2hyipmui2*mui[1:n1]+phyipmui)
  p2likb0b0Temp1 <- ( p2pipmui2*mui[(n1+1):n] + ppipmui )/(pii-1) - ppipmui^2*mui[(n1+1):n]/(1-pii)^2
  
  p2likb0b0Temp  <- c(p2likb0b0Temp0, p2likb0b0Temp1)
  
  for (i in 1:pdim) {
    for (j in i:pdim) {
      Jbar[i,j] <- t(p2likb0b0Temp)%*%(xij[,i]*xij[,j]*mui)
      if (i != j) Jbar[j,i] <- Jbar[i,j]
    }
  }
  
  if ( (sum(g0G0==0))==n1 ) {
    p2likg0b0Temp0 <- rep(0,n1)
  } else {
    p2likg0b0Temp0 <- (g0pG0-g0G0^2)*phyiptaui*phyipmui+g0G0*p2hyiptaumui
  }  
  p2likg0b0Temp1 <- p2piptaumui/(pii-1) - ppiptaui*ppipmui/(1-pii)^2
  p2likg0b0Temp <- c(p2likg0b0Temp0, p2likg0b0Temp1)
  for (i in 1:pdim) {
    for (j in (pdim+1):(pdim+qdim)) {
      Jbar[i,j] <- t(p2likg0b0Temp)%*%(xij[,i]*wij[,j-pdim]*mui)
      if (i != j) Jbar[j,i] <- Jbar[i,j]
    }
  }
  
  if ( (sum(g0G0==0))==n1 ) {
    p2likg0g0Temp0 <- rep(0,n1)
  } else {
    p2likg0g0Temp0 <- (g0pG0-g0G0^2)*phyiptaui^2+g0G0*p2hyiptaui2
  }
  
  p2likg0g0Temp1 <- -(p2piptaui2/(1-pii) + ppiptaui^2/(1-pii)^2)
  p2likg0g0Temp <- c(p2likg0g0Temp0, p2likg0g0Temp1)
  
  for (i in (pdim+1):(pdim+qdim)) {
    for (j in i:(pdim+qdim)) {
      Jbar[i,j] <- t(p2likg0g0Temp)%*%(wij[,i-pdim]*wij[,j-pdim])
      if (i != j) Jbar[j,i] <- Jbar[i,j]
    }
  } 
  
  Jbar <- -Jbar
  
  return(list(score, Jbar))
} # done
