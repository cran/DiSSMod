nr.nbinom <- function(xij, wij, theta0, yi, alpha0, select.dist, eps, itmax, trunc.num) {
  n <- nrow(xij); pdim <- ncol(xij); qdim <- ncol(wij)
  n1 <- n-sum(is.nan(yi))
  if(n<=n1) return(list(theta<-NA, Jbar<-NA, likval<-NA))
  pardim <- pdim+qdim+1      # +1 par.
  
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
  } else if (select.dist=="logistic"){
    fga <- function(ga02) {
      temp1 <- sum(plogis( ga02%*%t(wij[1:n1,]) ,log.p=TRUE) )
      temp2 <- sum(plogis( ga02%*%t(wij[(n1+1):n,]) , lower.tail=FALSE, log.p=TRUE ) )
      return(loglik <- -(temp1+temp2))
    }
  } # same
  theta <- theta0
  
  if (is.null(theta0)) # To give initial values
  {
    observed <- !is.na(yi)
    b0b1glm <- glm.nb(yi[1:n1] ~ xij[1:n1,2:pdim]) # nbinom_glm
    estpsi <- b0b1glm$theta
#    estpsi <- theta.md(yi[1:n1], fitted(b0b1glm), dfr = df.residual(b0b1glm))
#    estpsi <- theta.ml(yi[1:n1], fitted(b0b1glm))
#    estpsi <- theta.mm(yi[1:n1], fitted(b0b1glm), dfr = df.residual(b0b1glm))
    ga_esti <- optim(c(rep(1,qdim)),fga)$par
#      if(select.dist=="gumbel") {
#        ga_esti <- glm(as.numeric(observed) ~ wij[,2:qdim], family=binomial(logit))
#      } else {
#        ga_esti <- glm(as.numeric(observed) ~ wij[,2:qdim], family=binomial(probit))
#      }
#    theta <- as.numeric(c(b0b1glm$coef, ga_esti$coef, log(estpsi)))    # initial values
    theta <- as.numeric(c(b0b1glm$coef, ga_esti, log(estpsi)))    # initial values
  }
  
  alpha <- alpha0
  diff <- 1
  it <- 0; initmax <- 100 # To avoid an infinite loop
  likval <- like_nbinom(xij, wij, theta, yi, alpha, select.dist, trunc.num)
  if(is.na(likval)) {
    return(list(theta<-NA, Jbar<-NA, likval<-NA)) # Key Point where NA occurs #
  }
  
  while ((diff > eps) && (it < itmax))
  { 
    it <- it+1
    theta.old <- theta; likval.old <- likval
    
    # score and observed information matrix ------------------------------------------------------
    score_Jbar <- eval_score_obsinf_nbinom(xij, wij, theta, yi, alpha, select.dist, trunc.num) 
    
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
      
      likval <- like_nbinom(xij, wij, theta, yi, alpha, select.dist, trunc.num)
      
      if(is.na(likval)) likval <- -Inf
      
#cat("likval, diff",likval,diff,"\n")
#cat("tehta",theta,"\n")

      delta <- delta/2
    } # end of while((likval < likval.old) && (init < initmax))
    diff <- sum(abs(theta-theta.old)) 
#cat("within while",it,likval,"\n")
  } # end of while ((diff > eps) && (it < itmax))
  list(theta<-theta, Jbar<-Jbar, likval<-likval) # estimated parameters and Fisher observed information matrix
} # end of nr.nbinom


like_nbinom <- function(xij, wij, theta, yi, alpha, select.dist, trunc.num) {
  # log-likelihood function ---------------------------------------------
  pdim <- ncol(xij); qdim <- ncol(wij)
  n <- nrow(xij)
  pardim <- pdim+qdim+1
  n1 <- n-sum(is.nan(yi)) # count # of NaN's
  
  xibeta <- colSums(diag(theta[1:pdim])%*%t(xij))
  wigamma <- colSums(diag(theta[(pdim+1):(pdim+qdim)])%*%t(wij))
  
  b <- theta[pardim]
  mui <- exp(xibeta)
  mui[ which(mui >= .Machine$double.xmax) ] <- .Machine$double.xmax
  
  mu_i <- exp(-xibeta)
  mu_i[ which(mu_i >= .Machine$double.xmax) ] <- .Machine$double.xmax
  
  
  alpha <- alpha
  taui <- wigamma
  taui[ which( taui <= -.Machine$double.xmax ) ] <- -.Machine$double.xmax
  taui[ which( taui >= .Machine$double.xmax ) ] <- .Machine$double.xmax
  
  etai <- alpha*exp(-xibeta)
  etai[ which( etai <= -.Machine$double.xmax) ] <- -.Machine$double.xmax
  etai[ which( etai >= .Machine$double.xmax ) ] <- .Machine$double.xmax
  
  #  etai <- alpha/mui
  k <- 0:trunc.num
  
  # range of exp_b: 0 < exp(-b) < Inf
  expb <- exp(b)
  if(expb >= .Machine$double.xmax) {
    expb <- .Machine$double.xmax
  }
  
  if(select.dist=="gumbel") {
    hyi <- exp(taui+etai*yi) 
    hyi[ which(hyi <= .Machine$double.eps) ] <- .Machine$double.eps
    hyi[ which(hyi >= .Machine$double.xmax) ] <- .Machine$double.xmax
  } else {
    hyi <- taui+etai*yi
    hyi[ which(hyi <= -.Machine$double.xmax) ] <- -.Machine$double.xmax
    hyi[ which(hyi >= .Machine$double.xmax) ] <- .Machine$double.xmax
  }
  
  pii <- c()
  tempLik1 <- c()
  
  if(select.dist=="gumbel"){
    prob <- 1/( 1 + exp(xibeta[1:n1]-b) )
    prob[ which( prob <= .Machine$double.eps ) ] <- .Machine$double.eps
    prob[ which( prob >= 1-.Machine$double.eps ) ] <- 1 - .Machine$double.eps
    tempLik0 <- dnbinom(yi[1:n1], expb, prob, log=TRUE ) + pexp(hyi[1:n1], log.p=TRUE)
#tempLik0 <- dnbinom(yi[1:n1], expb, expb/(mui[1:n1]+expb), log=TRUE ) + log( 1 - exp( -hyi[1:n1] ) )
    for(i in (n1+1):n) {
      # range of pro: 0 < pro < 1
      pro <- 1/( 1 + exp(xibeta[i]-b) )
      if ( pro <= .Machine$double.eps ) {
        pro <- .Machine$double.eps
      } else if ( pro >= 1-.Machine$double.eps ) {
        pro <- 1-.Machine$double.eps
      }
      tempLik1[i-n1] <- log(  sum( exp( dnbinom(k, expb, pro, log=TRUE) ) * exp( -exp( taui[i]+etai[i]*k ) ) ) + .Machine$double.eps )
    }
  } else if(select.dist=="normal") {
    prob <- 1/( 1 + exp(xibeta[1:n1]-b) )
    prob[ which( prob <= .Machine$double.eps ) ] <- .Machine$double.eps
    prob[ which( prob >= 1-.Machine$double.eps ) ] <- 1 - .Machine$double.eps
    tempLik0 <- dnbinom(yi[1:n1], expb, prob, log=TRUE ) + pnorm( hyi[1:n1] , log.p=TRUE) 
    for(i in (n1+1):n) {
      # range of pro: 0 < pro < 1
      pro <- 1/( 1 + exp(xibeta[i]-b) )
      if ( pro <= .Machine$double.eps ) {
        pro <- .Machine$double.eps
      } else if ( pro >= 1-.Machine$double.eps ) {
        pro <- 1-.Machine$double.eps
      }
      pi <- sum( exp( dnbinom(k, expb, pro, log=TRUE) ) * exp( pnorm( taui[i]+etai[i]*k, log.p=TRUE) ) )
      tempLik1[i-n1] <- 1-pi  + .Machine$double.eps
      if(tempLik1[i-n1] < 0) {tempLik1[i-n1] <- .Machine$double.eps}
      tempLik1[i-n1] <- log( tempLik1[i-n1] )
    }
  } else if(select.dist=="logistic") {
    prob <- 1/( 1 + exp(xibeta[1:n1]-b) )
    prob[ which( prob <= .Machine$double.eps ) ] <- .Machine$double.eps
    prob[ which( prob >= 1-.Machine$double.eps ) ] <- 1 - .Machine$double.eps
    tempLik0 <- dnbinom(yi[1:n1], expb, prob, log=TRUE ) + plogis( hyi[1:n1] , log.p=TRUE)
    for(i in (n1+1):n) {
      pro <- 1/( 1 + exp(xibeta[i]-b) )
      if ( pro <= .Machine$double.eps ) {
        pro <- .Machine$double.eps
      } else if ( pro >= 1-.Machine$double.eps ) {
        pro <- 1-.Machine$double.eps
      }
      pi <- sum( exp(dnbinom(k, expb, pro, log=TRUE) ) * exp( plogis( taui[i]+etai[i]*k , log.p=TRUE) ) )
      tempLik1[i-n1] <-  1-pi  + .Machine$double.eps
      
      if(tempLik1[i-n1] < 0) {tempLik1[i-n1] <- .Machine$double.eps}
      tempLik1[i-n1] <- log( tempLik1[i-n1] )
    }
  }
  tempLik <-  c(tempLik0, tempLik1)  
  
  return(sum(tempLik))
}




eval_score_obsinf_nbinom <- function(xij, wij, theta, yi, alpha, select.dist, trunc.num) {
  # score function ---------------------------------------------------
  n <- nrow(xij); pdim <- ncol(xij); qdim <- ncol(wij)
  pardim <- pdim+qdim+1
  n1 <- n-sum(is.nan(yi)) # count # of NaN's
  
  score <- c()
  pii <- c()
  ppipmui <- c()
  ppiptaui <- c()
  ppippsi <- c()
  
  xibeta <- colSums(diag(theta[1:pdim])%*%t(xij))
  wigamma <- colSums(diag(theta[(pdim+1):(pdim+qdim)])%*%t(wij))

  b <- theta[pardim]
  mui <- exp(xibeta)
  mui[ which(mui >= .Machine$double.xmax) ] <- .Machine$double.xmax
  
  mu_i <- exp(-xibeta)
  mu_i[ which(mu_i >= .Machine$double.xmax) ] <- .Machine$double.xmax
  
  alpha <- alpha
  taui <- wigamma
  taui[ which( taui <= -.Machine$double.xmax ) ] <- -.Machine$double.xmax
  taui[ which( taui >= .Machine$double.xmax ) ] <- .Machine$double.xmax
  
  etai <- alpha*exp(-xibeta)
  etai[ which( etai <= -.Machine$double.xmax ) ] <- -.Machine$double.xmax
  etai[ which( etai >= .Machine$double.xmax ) ] <- .Machine$double.xmax
  #  etai <- alpha/mui
  k <- 0:trunc.num
  
  exp_b <- exp(-b)
  #  if ( exp_b <= .Machine$double.eps ) {
  #    exp_b  <- .Machine$double.eps
  #  } else if ( exp_b >= .Machine$double.xmax ) {
  #    exp_b <- .Machine$double.xmax
  #  }
  if(exp_b >= .Machine$double.xmax) {
    exp_b <- .Machine$double.xmax
  }
  
  # range of exp(b) : 0 < exp(b) < Inf
  expb <- exp(b)
  #  if ( expb <= .Machine$double.eps ) {
  #    expb <- .Machine$double.eps
  #  } else if ( expb >= .Machine$double.xmax ) {
  #    expb <- .Machine$double.xmax
  #  }
  if(expb >= .Machine$double.xmax) {
    expb <- .Machine$double.xmax
  }
  
  
  #range of exp(b+xibeta): 0 < exp(b+xibeta) < Inf
  expbx <- exp(b+xibeta)
  #  expbx [ which( expbx <= .Machine$double.eps ) ] <- .Machine$double.eps
  expbx [ which( expbx >= .Machine$double.xmax ) ] <- .Machine$double.xmax
  
  if(select.dist=="gumbel") {
    hyi <- exp(taui+etai*yi) 
    hyi[ which( hyi <= .Machine$double.eps ) ] <- .Machine$double.eps
    hyi[ which( hyi >= .Machine$double.xmax ) ] <- .Machine$double.xmax
    g0G0 <- exp( dexp(hyi[1:n1],log=TRUE) - pexp(hyi[1:n1], log.p=TRUE) )
    g0G0[ which( g0G0 >= .Machine$double.xmax ) ] <- .Machine$double.xmax
    g0pG0 <- -g0G0
    phyipmui <- -etai[1:n1]*mu_i[1:n1]*yi[1:n1]*hyi[1:n1]
    phyipmui[ which( phyipmui <= -.Machine$double.xmax ) ] <- -.Machine$double.xmax
    phyipmui[ which( phyipmui >= .Machine$double.xmax ) ] <- .Machine$double.xmax
    phyiptaui <- hyi[1:n1]
  } else if(select.dist=="normal") {
    hyi <- taui+etai*yi
    hyi[ which(hyi <= -.Machine$double.xmax) ] <- -.Machine$double.xmax
    hyi[ which(hyi >= .Machine$double.xmax) ] <- .Machine$double.xmax
    # First calculate log values and then exponentiate    
    g0G0 <- exp( dnorm(hyi[1:n1],log=TRUE) - pnorm(hyi[1:n1],log.p=TRUE) )
    g0G0[ which( g0G0 >= .Machine$double.xmax ) ] <- .Machine$double.xmax
    g0pG0 <- -hyi[1:n1]*g0G0
    g0pG0[ which( g0pG0 <= -.Machine$double.xmax ) ] <- -.Machine$double.xmax
    g0pG0[ which( g0pG0 >= .Machine$double.xmax ) ] <- .Machine$double.xmax
    phyipmui <- -etai[1:n1]*yi[1:n1]*mu_i[1:n1]
    phyiptaui <- rep(1,n1)
  } else if(select.dist=="logistic") {
    hyi <- taui+etai*yi
    hyi[ which(hyi <= -.Machine$double.xmax) ] <- -.Machine$double.xmax
    hyi[ which(hyi >= .Machine$double.xmax) ] <- .Machine$double.xmax
    phyipmui <- -etai[1:n1]*yi[1:n1]*mu_i[1:n1]
    phyiptaui <- rep(1,n1)
    g0G0 <- plogis(-(hyi[1:n1]))
    g0pG0 <- (2*plogis(-hyi[1:n1])-1)*g0G0
  }
  
  if (select.dist=="gumbel") {
    for(i in (n1+1):n) {
      prob <- 1/( 1 + exp(xibeta[i]-b) )
      prob[ which( prob <= .Machine$double.eps ) ] <- .Machine$double.eps
      prob[ which( prob >= 1-.Machine$double.eps ) ] <- 1 - .Machine$double.eps
      exptek <- exp(taui[i]+etai[i]*k)
      exptek[ which( exptek >= .Machine$double.xmax ) ] <- .Machine$double.xmax
      
      pii[i-n1] <- 1 - sum( exp( dnbinom(k, expb, prob, log=TRUE) ) * exp(-exp(taui[i]+etai[i]*k)) )
#pii[i-n1] <- 1 - sum( dnbinom(k, expb, expb/(mui[i]+expb)) * exp(-exp(taui[i]+etai[i]*k)) )
      ppipmui[i-n1] <- - sum( exp( dnbinom(k, expb, prob, log=TRUE) ) * exp(-exp(taui[i]+etai[i]*k)) *
                             ( -expb/(mui[i]+expb) + expb*k*mu_i[i]/(mui[i]+expb) + etai[i]*k*mu_i[i]*exptek ) )
      ppiptaui[i-n1] <- sum ( exp( dnbinom(k, expb, prob, log=TRUE) ) * exptek * exp(-exp(taui[i]+etai[i]*k)) )
      ppippsi[i-n1] <- - sum( exp( dnbinom(k, expb, prob, log=TRUE) ) * exp(-exp(taui[i]+etai[i]*k)) * expb * 
                             ( digamma(k+expb) - digamma(expb) + b - log(mui[i]+expb) + (mui[i]-k)/(mui[i]+expb) ) )
    }
  } else if(select.dist=="normal"){
    for (i in (n1+1):n) {
      prob <- 1/( 1 + exp(xibeta[i]-b) )
      prob[ which( prob <= .Machine$double.eps ) ] <- .Machine$double.eps
      prob[ which( prob >= 1-.Machine$double.eps ) ] <- 1 - .Machine$double.eps
      pii[i-n1] <- sum( exp( dnbinom(k,expb,prob, log=TRUE) + pnorm(taui[i]+etai[i]*k, log.p=TRUE) ) )
      ppipmui[i-n1] <- sum( exp( dnbinom(k,expb,prob, log=TRUE))* ( pnorm(taui[i]+etai[i]*k)*expb/(mui[i]+expb)*(k*mu_i[i]-1) - 
                                                                     dnorm(taui[i]+etai[i]*k)*etai[i]*k*mu_i[i] ) )
      ppiptaui[i-n1] <- sum( exp( dnbinom(k,expb,prob, log=TRUE))*dnorm(taui[i]+etai[i]*k) )
      ppippsi[i-n1] <- sum( exp( dnbinom(k,expb,prob, log=TRUE))*pnorm(taui[i]+etai[i]*k)*expb*
                           ( digamma(k+expb) - digamma(expb) + (mui[i]-k)/(mui[i]+expb)+b-log(mui[i]+expb) ) )
    }
  } else if(select.dist=="logistic") {
    for(i in (n1+1):n) {
      prob <- 1/( 1 + exp(xibeta[i]-b) )
      prob[ which( prob <= .Machine$double.eps ) ] <- .Machine$double.eps
      prob[ which( prob >= 1-.Machine$double.eps ) ] <- 1 - .Machine$double.eps
      pii[i-n1] <- sum( exp( dnbinom(k,expb,prob, log=TRUE) + plogis(taui[i]+etai[i]*k, log.p=TRUE) ) )
      ppipmui[i-n1] <- sum( exp( dnbinom(k,expb,prob, log=TRUE) ) * plogis(taui[i]+etai[i]*k) *
                              ( prob*(k*mu_i[i]-1) - etai[i]*k*mu_i[i]*plogis(-(taui[i]+etai[i]*k)) ) )
      ppiptaui[i-n1] <- sum( exp( dnbinom(k, expb, prob, log=TRUE) ) * plogis(taui[i]+etai[i]*k) * plogis(-(taui[i]+etai[i]*k)) ) 
      ppippsi[i-n1] <- sum( exp( dnbinom(k, expb, prob, log=TRUE) ) * plogis(taui[i]+etai[i]*k) * expb *
                              ( digamma(k+expb) - digamma(expb) + b - log(mui[i]+expb) + (mui[i]-k)/(mui[i]+expb) ) )
    }
  }
  
  plfpmui <- yi*mu_i - (expb+yi)/(mui+expb)
  plfppsi <- expb * ( digamma(yi+expb) - digamma(expb) + b - log(mui+expb) + (mui-yi)/(mui+expb) )
  
  pii[which(pii >= 1-.Machine$double.eps)] <- 1-.Machine$double.eps
  di0pii <- -1/(1-pii) 
  di0pii2 <- - 1/(1-pii)^2 
  
  
  for (j in 1:pdim) {
    score[j] <- t( (plfpmui[1:n1] + g0G0*phyipmui)*mui[1:n1] )%*%xij[1:n1,j] +
      t( di0pii*ppipmui* mui[(n1+1):n] )%*%xij[(n1+1):n,j]
  }
  
  for (j in 1:qdim) {
    if ( (sum(g0G0==0))==n1 ) {
      score[pdim+j] <- t(di0pii*ppiptaui)%*%wij[(n1+1):n,j] 
    } else {
      score[pdim+j] <- t(g0G0*phyiptaui)%*%wij[1:n1,j] + t(di0pii*ppiptaui)%*%wij[(n1+1):n,j] 
    }  
  }

  score[pdim+qdim+1] <- sum(plfppsi[1:n1]) + sum(di0pii*ppippsi)
  
  
  # observed information matrix  ------------------------------------------------------
  Jbar <- matrix(0,pardim,pardim)

  p2pipmui2 <- c() 
  p2piptaui2 <- c()
  p2pippsi2 <- c()
  p2piptaumui <- c()
  p2piptaupsi <- c()
  p2pippsimui <- c()
  
  
  if (select.dist=="gumbel") {
    p2hyipmui2 <- etai[1:n1]*yi[1:n1]*mu_i[1:n1]^2*(2+etai[1:n1]*yi[1:n1])*hyi[1:n1]
    p2hyiptaumui <- -etai[1:n1]*yi[1:n1]*mu_i[1:n1]*hyi[1:n1]
    p2hyiptaui2 <- hyi[1:n1] 
  } else if(select.dist=="normal" | select.dist=="logistic"){
    p2hyipmui2 <- 2*etai[1:n1]*yi[1:n1]*mu_i[1:n1]^2
    p2hyiptaumui <- rep(0,n1)
    p2hyiptaui2 <- rep(0,n1)
  }
  
  
  if (select.dist=="gumbel") {
    for(i in (n1+1):n) {
      prob <- 1/( 1 + exp(xibeta[i]-b) )
      prob[ which( prob <= .Machine$double.eps ) ] <- .Machine$double.eps
      prob[ which( prob >= 1-.Machine$double.eps ) ] <- 1 - .Machine$double.eps
      exptek <- exp(taui[i]+etai[i]*k)
      exptek[ which( exptek >= .Machine$double.xmax ) ] <- .Machine$double.xmax
      tem <- ( expb*(k*mu_i[i]-1)/(mui[i]+expb) + etai[i]*k*mu_i[i]*exptek )^2
      tem[ which( tem >= .Machine$double.xmax ) ] <- .Machine$double.xmax
      
      p2pipmui2[i-n1] <- - sum( exp( dnbinom(k, expb, prob, log=TRUE) ) * exp(-exp(taui[i]+etai[i]*k)) * 
                               ( tem + expb*( 1 - k*(2*mui[i]+expb)*mu_i[i]^2 )/(mui[i]+expb)^2 -
                                   etai[i]*k*exptek*(2+etai[i]*k)*mu_i[i]^2 ) )
      p2piptaumui[i-n1] <- sum( exp( dnbinom(k, expb, prob, log=TRUE) ) * exp(-exp(taui[i]+etai[i]*k)) * exptek *
                               ( expb*(k*mu_i[i]-1)/(mui[i]+expb) + etai[i]*k*mu_i[i]*(exptek-1) ) )
      p2pippsimui[i-n1] <- -sum( exp( dnbinom(k, expb, prob, log=TRUE) ) * exp(-exp(taui[i]+etai[i]*k)) * expb *
                                ( (expb*(k*mu_i[i]-1)/(mui[i]+expb) + etai[i]*k*mu_i[i]*exptek) *
                                    (  digamma(k+expb) - digamma(expb) + b - log(mui[i]+expb) + (mui[i]-k)/(mui[i]+expb) ) +
                                    (k-mui[i])/(mui[i]+expb)^2 ) )
      p2piptaui2[i-n1] <- sum( exp( dnbinom(k, expb, prob, log=TRUE) ) * exp(-exp(taui[i]+etai[i]*k)) * exptek *
                              ( 1-exptek ) )
      p2piptaupsi[i-n1] <- sum( exp( dnbinom(k, expb, prob, log=TRUE) ) * exp(-exp(taui[i]+etai[i]*k)) * exptek * expb *
                               ( digamma(k+expb) - digamma(expb) + b - log(mui[i]+expb) + (mui[i]-k)/(mui[i]+expb) ) )
      p2pippsi2[i-n1] <- -sum( exp( dnbinom(k, expb, prob, log=TRUE) ) * exp(-exp(taui[i]+etai[i]*k)) * 
                              ( expb * (digamma(k+expb) - digamma(expb) + b - log(mui[i]+expb) + (mui[i]-k)/(mui[i]+expb) ) +
                                  exp(2*b) * (digamma(k+expb) - digamma(expb) + b - log(mui[i]+expb) + (mui[i]-k)/(mui[i]+expb))^2 +
                                  exp(2*b)*( trigamma(k+expb) - trigamma(expb) + mui[i]*exp(-b)/(mui[i]+expb) - (mui[i]-k)/(mui[i]+expb)^2 ) ) )
      }
  } else if(select.dist=="normal") {
    for(i in (n1+1):n) {
      prob <- 1/( 1 + exp(xibeta[i]-b) )
      prob[ which( prob <= .Machine$double.eps ) ] <- .Machine$double.eps
      prob[ which( prob >= 1-.Machine$double.eps ) ] <- 1 - .Machine$double.eps
      p2pipmui2[i-n1] <- sum( exp( dnbinom(k,expb,prob, log=TRUE))*( pnorm(taui[i]+etai[i]*k)*expb/(mui[i]+expb)^2*( expb*(k*mu_i[i]-1)^2+
                                                                                                                          1-k*(2*mui[i]+expb)*mu_i[i]^2 )-
                                                                      dnorm(taui[i]+etai[i]*k)*etai[i]*k*mu_i[i]*( expb*2*(k*mu_i[i]-1)/(mui[i]+expb)-
                                                                                                                     2*mu_i[i]+etai[i]*k*mu_i[i]*(taui[i]+etai[i]*k) ) ) )
      p2piptaumui[i-n1] <- sum( exp( dnbinom(k,expb,prob, log=TRUE))*dnorm(taui[i]+etai[i]*k)*
                               ( expb*(k*mu_i[i]-1)/(mui[i]+expb)+etai[i]*k*mu_i[i]*(taui[i]+etai[i]*k) ) ) 
      p2pippsimui[i-n1] <- sum( exp( dnbinom(k,expb,prob, log=TRUE))* ( pnorm(taui[i]+etai[i]*k)*expb*(k*mu_i[i]-1)/(mui[i]+expb)*
                                                                         ( expb*( digamma(k+expb)-digamma(expb)
                                                                                    +(mui[i]-k)/(mui[i]+expb)+b-log(mui[i]+expb) ) +
                                                                             mui[i]/(mui[i]+expb) ) -
                                                                         etai[i]*k*mu_i[i]*dnorm(taui[i]+etai[i]*k)*expb*
                                                                         ( digamma(k+expb)-digamma(expb)+(mui[i]-k)/(mui[i]+expb)+
                                                                             b-log(mui[i]+expb) ) ) )
      p2piptaui2[i-n1] <- -sum( exp( dnbinom(k,expb,prob, log=TRUE))*(taui[i]+etai[i]*k)*dnorm(taui[i]+etai[i]*k) )
      p2piptaupsi[i-n1] <- sum( exp( dnbinom(k,expb,prob, log=TRUE))*dnorm(taui[i]+etai[i]*k)*expb*
                               ( digamma(k+expb)-digamma(expb)+(mui[i]-k)/(mui[i]+expb)+b-log(mui[i]+expb) ) )
      p2pippsi2[i-n1] <- sum( exp( dnbinom(k,expb,prob, log=TRUE))*pnorm(taui[i]+etai[i]*k)*expb*
                             ( expb*( digamma(k+expb)-digamma(expb)+(mui[i]-k)/(mui[i]+expb)+b-log(mui[i]+expb) )^2+
                                 ( digamma(k+expb)-digamma(expb)+(mui[i]-k)/(mui[i]+expb)+b-log(mui[i]+expb) )+
                                 expb*( trigamma(k+expb)-trigamma(expb)+(k-mui[i])/(mui[i]+expb)^2+mui[i]*exp(-b)/(mui[i]+expb) ) ) )
    }
  } else if(select.dist=="logistic") {
    for(i in (n1+1):n) {
      prob <- 1/( 1 + exp(xibeta[i]-b) )
      prob[ which( prob <= .Machine$double.eps ) ] <- .Machine$double.eps
      prob[ which( prob >= 1-.Machine$double.eps ) ] <- 1 - .Machine$double.eps
      p2pipmui2[i-n1] <- sum( exp( dnbinom(k, expb, prob, log=TRUE) ) * plogis( taui[i]+etai[i]*k) * 
                             ( ( expb/(mui[i]+expb)*(k*mu_i[i]-1) - etai[i]*k*mu_i[i]*plogis(-(taui[i]+etai[i]*k)) )^2 +
                                 expb/(expb+mui[i])^2*(1-k*(2*mui[i]+expb)*mu_i[i]^2) +
                                 etai[i]*k*mu_i[i]^2*plogis(-(taui[i]+etai[i]*k))*(2-etai[i]*k*plogis(taui[i]+etai[i]*k))))
      p2piptaumui[i-n1] <- sum( exp( dnbinom(k,expb,prob, log=TRUE) ) * plogis(taui[i]+etai[i]*k) * 
                               ( plogis(-(taui[i]+etai[i]*k)) * ( expb/(mui[i]+expb)*(k*mu_i[i]-1) - etai[i]*k*mu_i[i]*
                                                                    (plogis(-(taui[i]+etai[i]*k)) - plogis(taui[i]+etai[i]*k)) ) ) )
      p2pippsimui[i-n1] <- sum( exp( dnbinom(k,expb,prob, log=TRUE) ) * plogis(taui[i]+etai[i]*k) * expb * 
                               ( ( expb/(mui[i]+expb)*(k*mu_i[i]-1) - etai[i]*k*mu_i[i]*plogis(-(taui[i]+etai[i]*k)) ) *
                                   ( digamma(k+expb) - digamma(expb) + b - log(mui[i]+expb) + (mui[i]-k)/(mui[i]+expb) ) +
                                   mui[i]/(mui[i]+expb)^2*(k*mu_i[i]-1) ) )
      p2piptaui2[i-n1] <- sum( exp( dnbinom(k,expb,prob, log=TRUE) ) * plogis(taui[i]+etai[i]*k) * plogis(-(taui[i]+etai[i]*k)) *
                              ( plogis(-(taui[i]+etai[i]*k)) - plogis(taui[i]+etai[i]*k) ) )
      p2piptaupsi[i-n1] <- sum( exp( dnbinom(k,expb,prob, log=TRUE) ) * plogis(taui[i]+etai[i]*k) * expb *
                               plogis(-(taui[i]+etai[i]*k)) * ( digamma(k+expb) - digamma(expb) + b - log(mui[i]+expb) +
                                                                  (mui[i]-k)/(mui[i]+expb) ) ) 
      p2pippsi2[i-n1] <- sum( exp( dnbinom(k,expb,prob, log=TRUE) ) * plogis(taui[i]+etai[i]*k) *
                             ( expb^2 * ( digamma(k+expb)-digamma(expb)+b-log(mui[i]+expb)+(mui[i]-k)/(mui[i]+expb) )^2 +
                                 expb * ( digamma(k+expb)-digamma(expb)+b-log(mui[i]+expb)+(mui[i]-k)/(mui[i]+expb) ) +
                                 expb^2 * ( trigamma(k+expb)-trigamma(expb)+mui[i]/(mui[i]+expb)/expb-(mui[i]-k)/(mui[i]+expb)^2 ) ) )
      
    }
  }
  
  p2lfpmui2 <- -yi*mu_i^2 + (expb+yi)/(mui+expb)^2
  p2lfppsi2 <- expb*( digamma(expb+yi) - digamma(expb) + b - log(mui+expb) + (mui-yi)/(mui+expb) ) +
    exp(2*b)*( trigamma(yi+expb) - trigamma(expb) + mui*exp(-b)/(mui+expb) - (mui-yi)/(mui+expb)^2 )
  p2lfppsimui <- expb*(yi-mui)/(mui+expb)^2
  
  
  p2likb0b0Temp0 <- (p2lfpmui2[1:n1]+g0pG0*phyipmui*phyipmui-g0G0^2*phyipmui*phyipmui+g0G0*p2hyipmui2)*mui[1:n1] + (plfpmui[1:n1]+g0G0*phyipmui)
#p2likb0b0Temp0 <- (p2lfpmui2[1:n1]+(g0pG0-g0G0^2)*phyipmui^2+g0G0*p2hyipmui2)*mui[1:n1]^2 + (plfpmui[1:n1]+g0G0*phyipmui)*mui[1:n1] 
  p2likb0b0Temp1 <-  di0pii2*( p2pipmui2*(1-pii) + ppipmui*ppipmui )*mui[(n1+1):n]+ di0pii*ppipmui
  
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
    p2likg0b0Temp0 <- g0pG0*phyiptaui*phyipmui-g0G0^2*phyiptaui*phyipmui + g0G0*p2hyiptaumui
#p2likg0b0Temp0 <- (g0pG0-g0G0^2)*phyiptaui*phyipmui + g0G0*p2hyiptaumui
  }  
  p2likg0b0Temp1 <- di0pii2*( p2piptaumui*(1-pii)+ppiptaui*ppipmui )
  
  p2likg0b0Temp <- c(p2likg0b0Temp0, p2likg0b0Temp1)
  
  for (i in 1:pdim) {
    for (j in (pdim+1):(pdim+qdim)) {
      Jbar[i,j] <- t(p2likg0b0Temp)%*%(xij[,i]*wij[,j-pdim]*mui)
      if (i != j) Jbar[j,i] <- Jbar[i,j]
    }
  }
  
  
  p2likp0b0Temp0 <- p2lfppsimui[1:n1]
  p2likp0b0Temp1 <-di0pii2*( p2pippsimui*(1-pii) + ppipmui*ppippsi )
  
  p2likp0b0Temp <- c(p2likp0b0Temp0, p2likp0b0Temp1)
  
  for (i in 1:pdim) {
    j <- pardim
    Jbar[i,j] <- t(p2likp0b0Temp)%*%(xij[,i]*mui)
    if(i != j) Jbar[j,i] <- Jbar[i,j]
  }
  
  
  if ( (sum(g0G0==0))==n1 ) {
    p2likg0g0Temp0 <- rep(0,n1)
  } else {
    p2likg0g0Temp0 <- g0pG0*phyiptaui*phyiptaui-g0G0^2*phyiptaui*phyiptaui + g0G0*p2hyiptaui2
#p2likg0g0Temp0 <- (g0pG0-g0G0^2)*phyiptaui^2 + g0G0*p2hyiptaui2
  }
  p2likg0g0Temp1 <- di0pii2*( p2piptaui2*(1-pii) + ppiptaui^2 )
  
  p2likg0g0Temp <- c(p2likg0g0Temp0, p2likg0g0Temp1)
  
  for (i in (pdim+1):(pdim+qdim)) {
    for (j in i:(pdim+qdim)) {
      Jbar[i,j] <- t(p2likg0g0Temp)%*%(wij[,i-pdim]*wij[,j-pdim])
      if (i != j) Jbar[j,i] <- Jbar[i,j]
    }
  } 

  
  p2likp0g0Temp0 <- rep(0,n1) 
  p2likp0g0Temp1 <-  di0pii2*( p2piptaupsi*(1-pii) + ppiptaui*ppippsi )
  
  p2likp0g0Temp <- c(p2likp0g0Temp0, p2likp0g0Temp1)
  
  for (i in (pdim+1):(pdim+qdim)) {
    j <- pardim
    Jbar[i,j] <- t(p2likp0g0Temp)%*%wij[,i-pdim]
    if(i != j) Jbar[j,i] <- Jbar[i,j]
  }
  
  
  p2likp0p0Temp0 <- p2lfppsi2[1:n1]
  p2likp0p0Temp1 <- di0pii2*(p2pippsi2*(1-pii) + ppippsi^2)
  
  p2likp0p0Temp <- c(p2likp0p0Temp0, p2likp0p0Temp1)
  
  i <- j <- pardim
  Jbar[i,j] <- sum(p2likp0p0Temp)
  
  
  
  Jbar <- -Jbar
  
  return(list(score, Jbar))
}

