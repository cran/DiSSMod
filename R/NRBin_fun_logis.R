#--------------------------------------------------------------
#
# Newton Raphson Function for Binomial case (T ~ logistic)
#
# Programmed by Sang Kyu Lee at 12-19-2018
#
#--------------------------------------------------------------
eval_score_obsinf_bin <- function(mi, xij, wij, theta, yi, alpha, select.dist) {
  # score function ---------------------------------------------------(ok)
  n <- nrow(xij); pdim <- ncol(xij); qdim <- ncol(wij)
  pardim <- pdim+qdim
  n1 <- n-sum(is.na(yi)) # count # of Na's

  score <- c()

  xibeta <- colSums(diag(theta[1:pdim])%*%t(xij))
  wigamma <- colSums(diag(theta[(pdim+1):(pdim+qdim)])%*%t(wij))

  probi <- 1/(1+exp(-xibeta))
  probi[ which(probi < .Machine$double.eps) ] <- .Machine$double.eps
  probi[ which(probi > 1) ] <- 1 - .Machine$double.eps

  mui <- mi*probi

  alpha <- alpha
  taui <- wigamma
  taui[ which( taui <= -.Machine$double.xmax ) ] <- -.Machine$double.xmax
  taui[ which( taui >= .Machine$double.xmax ) ] <- .Machine$double.xmax

  etai <- alpha/mui
  etai[ which( etai <= -.Machine$double.xmax) ] <- -.Machine$double.xmax
  etai[ which( etai >= .Machine$double.xmax ) ] <- .Machine$double.xmax


  if(select.dist=="gumbel") {
    hyi <- exp(taui+etai*yi)
    g0G0 <- exp( dexp(hyi[1:n1],log=TRUE) - pexp(hyi[1:n1], log.p=TRUE) )
    #    g0G0 <- 1/(exp(hyi[1:n1]) - 1)
    #exp(-hyi[1:n1])/(1-exp(-hyi[1:n1]))
    g0pG0 <- -g0G0
  } else if(select.dist=="normal") {
    hyi <- taui+etai*yi
    #    g0G0 <- dnorm(hyi[1:n1])/pnorm(hyi[1:n1])
    # First calculate log values and then exponentiate
    g0G0 <- exp( dnorm(hyi[1:n1],log=TRUE)-pnorm(hyi[1:n1],log.p=TRUE) )
    g0pG0 <- -hyi[1:n1]*g0G0
  } else if(select.dist=="logistic") { # new
    hyi <- taui+etai*yi
    hyi[ which(hyi <= -.Machine$double.xmax) ] <- -.Machine$double.xmax
    hyi[ which(hyi >= .Machine$double.xmax) ] <- .Machine$double.xmax

    g0G0 <- exp( dlogis(hyi[1:n1],log=TRUE)-plogis(hyi[1:n1],log.p=TRUE) )
    g0pG0 <- -exp(-hyi[1:n1])*plogis(hyi[1:n1]) + 2*exp(-2*hyi[1:n1])*(plogis(hyi[1:n1])^2)
  }

  if (select.dist=="gumbel") {
    phyipmui <- -etai[1:n1]*yi[1:n1]/mui[1:n1]*hyi[1:n1]
    pii <- 1-(1-probi[(n1+1):n])*exp(-exp(taui[(n1+1):n]))-probi[(n1+1):n]*exp(-exp(taui[(n1+1):n]+alpha/probi[(n1+1):n]))
    ppipmui <- exp(-exp(taui[(n1+1):n]))-exp(-exp(taui[(n1+1):n]+alpha/probi[(n1+1):n]))*(1+alpha*exp(taui[(n1+1):n]+alpha/probi[(n1+1):n])/probi[(n1+1):n])
  } else if (select.dist=="normal") {
    phyipmui <- -etai[1:n1]*yi[1:n1]/mui[1:n1]
    pii <- (1-probi[(n1+1):n])*pnorm(taui[(n1+1):n])+probi[(n1+1):n]*pnorm(taui[(n1+1):n]+alpha/probi[(n1+1):n])
    ppipmui <- -pnorm(taui[(n1+1):n])+pnorm(taui[(n1+1):n]+alpha/probi[(n1+1):n])-dnorm(taui[(n1+1):n]+alpha/probi[(n1+1):n])*alpha/probi[(n1+1):n]
  } else if (select.dist=="logistic") { # new
    phyipmui <- -etai[1:n1]*yi[1:n1]/mui[1:n1]
    pii <- (1-probi[(n1+1):n])/(1+exp(-taui[(n1+1):n])) + probi[(n1+1):n]/(1+exp(-(taui[(n1+1):n] + alpha/probi[(n1+1):n])))
    ppipmui <- -1/(1 + exp(-taui[(n1+1):n])) + 1/(1+ exp(-(taui[(n1+1):n] + alpha/probi[(n1+1):n]))) -
      alpha * exp(-(taui[(n1+1):n] + alpha/probi[(n1+1):n])) / ( probi[(n1+1):n] * ( 1 + exp(-(taui[(n1+1):n] + alpha/probi[(n1+1):n])) )^2 )
  }

  for (j in 1:pdim) {
    score[j] <- t( ((yi[1:n1]-mui[1:n1])/(mui[1:n1]*(1-probi[1:n1])) + g0G0*phyipmui)*( mui[1:n1]*(mi[1:n1]-mui[1:n1])/mi[1:n1]))%*%xij[1:n1,j] +
      t( (ppipmui/(pii-1))*( mui[(n1+1):n]*(mi[(n1+1):n]-mui[(n1+1):n])/mi[(n1+1):n]) )%*%xij[(n1+1):n,j]
  }

  if (select.dist=="gumbel") {
    phyiptaui <- hyi[1:n1]
    ppiptaui <- (1-probi[(n1+1):n])*exp(taui[(n1+1):n])*exp(-exp(taui[(n1+1):n]))
    ppiptaui <- ppiptaui + probi[(n1+1):n]*exp(taui[(n1+1):n]+alpha/probi[(n1+1):n])*exp(-exp(taui[(n1+1):n]+alpha/probi[(n1+1):n]))
  } else if (select.dist=="normal") {
    phyiptaui <- rep(1,n1)
    ppiptaui <- (1-probi[(n1+1):n])*dnorm(taui[(n1+1):n])+probi[(n1+1):n]*dnorm(taui[(n1+1):n]+alpha/probi[(n1+1):n])
  } else if (select.dist=="logistic") { # new
    phyiptaui <- rep(1,n1)
    ppiptaui <- (1-probi[(n1+1):n])*exp(-taui[(n1+1):n])/( (1+exp(-taui[(n1+1):n]))^2 ) +
      probi[(n1+1):n]*exp(-(taui[(n1+1):n] + alpha/probi[(n1+1):n]))/( (1+exp(-(taui[(n1+1):n] + alpha/probi[(n1+1):n])))^2 )
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

  if (select.dist=="gumbel") {
    p2hyipmui2 <- etai[1:n1]*yi[1:n1]/mui[1:n1]^2*hyi[1:n1]*(2+etai[1:n1]*yi[1:n1])
    p2pipmui2 <- alpha^2/probi[(n1+1):n]^3*exp(taui[(n1+1):n]+alpha/probi[(n1+1):n])
    p2pipmui2 <- p2pipmui2*exp(-exp(taui[(n1+1):n]+alpha/probi[(n1+1):n]))*(1-exp(taui[(n1+1):n]+alpha/probi[(n1+1):n]))
  } else if (select.dist=="normal") {
    p2hyipmui2 <- 2*etai[1:n1]*yi[1:n1]/mui[1:n1]^2
    p2pipmui2 <- -alpha^2/probi[(n1+1):n]^3 * (taui[(n1+1):n]+alpha/probi[(n1+1):n])*dnorm(taui[(n1+1):n]+alpha/probi[(n1+1):n])
  } else if (select.dist=="logistic") { # new
    p2hyipmui2 <- 2*etai[1:n1]*yi[1:n1]/((mui[1:n1])^2)
    p2pipmui2 <- 2*(alpha^2)/(probi[(n1+1):n]^3) * exp(-(2*(taui[(n1+1):n] + alpha/probi[(n1+1):n]))) *
     plogis(taui[(n1+1):n] + alpha/probi[(n1+1):n])^3 - (alpha^2)/(probi[(n1+1):n]^3) *
     exp(-(taui[(n1+1):n] + alpha/probi[(n1+1):n])) * plogis(taui[(n1+1):n] + alpha/probi[(n1+1):n])^2
    # star <- 1-alpha*exp(-(taui[(n1+1):n] + alpha/probi[(n1+1):n]))*plogis(taui[(n1+1):n] + alpha/probi[(n1+1):n])/probi[(n1+1):n]
    # p2pipmui2 <- -alpha*exp(-(taui[(n1+1):n]+alpha/probi[(n1+1):n]))*(plogis(taui[(n1+1):n] + alpha/probi[(n1+1):n])^2)/(probi[(n1+1):n]^2)*star +
    #   alpha*exp(-(taui[(n1+1):n]+alpha/probi[(n1+1):n]))*(plogis(taui[(n1+1):n] + alpha/probi[(n1+1):n])^2)/probi[(n1+1):n]^2 -
    #   alpha^2*exp(-(taui[(n1+1):n]+alpha/probi[(n1+1):n]))*(plogis(taui[(n1+1):n] + alpha/probi[(n1+1):n])^2)/probi[(n1+1):n]^3 +
    #   alpha^2*exp(-2*(taui[(n1+1):n]+alpha/probi[(n1+1):n]))*(plogis(taui[(n1+1):n] + alpha/probi[(n1+1):n])^3)/probi[(n1+1):n]^3 # same result
  }

  p2likb0b0Temp0 <- -1+(g0pG0-g0G0^2)*phyipmui^2*mui[1:n1]*(1-probi[1:n1])+ g0G0*( p2hyipmui2*mui[1:n1]*(1-probi[1:n1])+phyipmui*(1-2*probi[1:n1]) )
  ##  p2likb0b0Temp0 <- p2likb0b0Temp0 -( (yi[1:n1]-mui[1:n1])/(mui[1:n1]*(1-probi[1:n1])) + g0G0*phyipmui )*( 1-2*probi[1:n1]+mui[1:n1]*(1-probi[1:n1])*(-1/mui[1:n1]^2+1/(mi[1:n1]-mui[1:n1])^2)/(1/mui[1:n1]+1/(mi[1:n1]-mui[1:n1])) )

  p2likb0b0Temp1 <- ( p2pipmui2*mui[(n1+1):n]*(1-probi[(n1+1):n]) + ppipmui*(1-2*probi[(n1+1):n]) )/(pii-1) - (ppipmui^2*mui[(n1+1):n]*(1-probi[(n1+1):n]))/(1-pii)^2
  ##  p2likb0b0Temp1 <- p2likb0b0Temp1 + ( ppipmui/(1-pii) )*( 1-2*probi[(n1+1):n]+mui[(n1+1):n]*(1-probi[(n1+1):n])*(-1/mui[(n1+1):n]^2+1/(mi[(n1+1):n]-mui[(n1+1):n])^2)/(1/mui[(n1+1):n]+1/(mi[(n1+1):n]-mui[(n1+1):n])) )

  p2likb0b0Temp  <- c(p2likb0b0Temp0, p2likb0b0Temp1)

  for (i in 1:pdim) {
    for (j in i:pdim) {
      Jbar[i,j] <- t(p2likb0b0Temp)%*%(xij[,i]*xij[,j]/( mui*(1-probi) *( 1/mui + 1/(mi-mui))^2 ) )
      if (i != j) Jbar[j,i] <- Jbar[i,j]
    }
  }

  if (select.dist=="gumbel") {
    p2hyiptaumui <- -etai[1:n1]*yi[1:n1]*hyi[1:n1]/mui[1:n1]
    p2piptaumui <- -exp(taui[(n1+1):n])*exp(-exp(taui[(n1+1):n]))+exp(taui[(n1+1):n]+alpha/probi[(n1+1):n])*exp(-exp(taui[(n1+1):n]+alpha/probi[(n1+1):n]))*(1+alpha/probi[(n1+1):n]*(exp(taui[(n1+1):n]+alpha/probi[(n1+1):n])-1))
  } else if (select.dist=="normal") {
    p2hyiptaumui <- rep(0,n1)
    p2piptaumui <- -dnorm(taui[(n1+1):n])+dnorm(taui[(n1+1):n]+alpha/probi[(n1+1):n])*(1+alpha*(taui[(n1+1):n]+alpha/probi[(n1+1):n])/probi[(n1+1):n])
  } else if (select.dist=="logistic") { # new
    p2hyiptaumui <- rep(0,n1)
    p2piptaumui <- -exp(-taui[(n1+1):n])*(plogis(taui[(n1+1):n])^2) + exp(-(taui[(n1+1):n] + alpha/probi[(n1+1):n]))*
      (plogis(taui[(n1+1):n] + alpha/probi[(n1+1):n])^2) - 2*alpha/probi[(n1+1):n] * exp(-2*(taui[(n1+1):n] + alpha/probi[(n1+1):n])) *
      (plogis(taui[(n1+1):n] + alpha/probi[(n1+1):n])^3) + alpha/probi[(n1+1):n] * exp(-(taui[(n1+1):n] + alpha/probi[(n1+1):n])) *
      (plogis(taui[(n1+1):n] + alpha/probi[(n1+1):n])^2)
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
      Jbar[i,j] <- t(p2likg0b0Temp)%*%(xij[,i]*wij[,j-pdim]/(1/mui+1/(mi-mui)))
      if (i != j) Jbar[j,i] <- Jbar[i,j]
    }
  }

  if (select.dist=="gumbel") {
    p2hyiptaui2 <- hyi[1:n1]
    p2piptaui2 <- (1-probi[(n1+1):n])*exp(taui[(n1+1):n])*exp(-exp(taui[(n1+1):n]))*(1-exp(taui[(n1+1):n]))
    p2piptaui2 <- p2piptaui2 + probi[(n1+1):n]*exp(taui[(n1+1):n]+alpha/probi[(n1+1):n])*exp(-exp(taui[(n1+1):n]+alpha/probi[(n1+1):n]))*(1-exp(taui[(n1+1):n]+alpha/probi[(n1+1):n]))
  } else if (select.dist=="normal") {
    p2hyiptaui2 <- rep(0,n1)
    p2piptaui2 <- (probi[(n1+1):n]-1)*taui[(n1+1):n]*dnorm(taui[(n1+1):n])-probi[(n1+1):n]*(taui[(n1+1):n]+alpha/probi[(n1+1):n])*dnorm(taui[(n1+1):n]+alpha/probi[(n1+1):n])
  } else if (select.dist=="logistic") {
    p2hyiptaui2 <- rep(0,n1)
    p2piptaui2 <- -(1-probi[(n1+1):n]) * exp(-taui[(n1+1):n]) * (plogis(taui[(n1+1):n])^2) +
      2 * (1-probi[(n1+1):n]) * exp(-2*taui[(n1+1):n]) * (plogis(taui[(n1+1):n])^3) -
      probi[(n1+1):n] * exp(-(taui[(n1+1):n] + alpha/probi[(n1+1):n])) * ( plogis(taui[(n1+1):n] + alpha/probi[(n1+1):n])^2 ) +
      2 * probi[(n1+1):n] * exp(-2*(taui[(n1+1):n] + alpha/probi[(n1+1):n])) * (plogis(taui[(n1+1):n] + alpha/probi[(n1+1):n])^3)
  }

  if ( (sum(g0G0==0))==n1 ) {
    p2likg0g0Temp0 <- rep(0,n1)
  } else {
    p2likg0g0Temp0 <- (g0pG0-g0G0^2)*phyiptaui^2 + g0G0*p2hyiptaui2

  }

  p2likg0g0Temp1 <- - (p2piptaui2/(1-pii) + ppiptaui^2/(1-pii)^2)
  p2likg0g0Temp <- c(p2likg0g0Temp0, p2likg0g0Temp1)

  for (i in (pdim+1):(pdim+qdim)) {
    for (j in i:(pdim+qdim)) {
      Jbar[i,j] <- t(p2likg0g0Temp)%*%(wij[,i-pdim]*wij[,j-pdim])
      if (i != j) Jbar[j,i] <- Jbar[i,j]
    }
  }

  Jbar <- -Jbar

  return(list(score, Jbar))
}

#----------------------------------------------------------------------------------
nr.bin <- function(mi, xij, wij, theta0, yi, alpha0, select.dist, eps, itmax) {
  n <- nrow(xij); pdim <- ncol(xij); qdim <- ncol(wij)
  n1 <- n-sum(is.nan(yi))
  if(n<=n1) return(list(theta<-NA, Jbar<-NA, likval<-NA))
  pardim <- pdim+qdim

  if(select.dist=="gumbel") {
    fga <- function(ga02) {
      temp1 <- sum(pexp( exp(ga02%*%t(wij[1:n1,])) , log.p=TRUE))
      #log(1-exp(-exp( ga02%*%t(wij[1:n1,]) ))))
      temp2 <- -sum(exp(  ga02%*%t(wij[(n1+1):n,]) ))
      return(loglik <- -(temp1+temp2))
    }
  } else if(select.dist=="normal") {
    fga <- function(ga02) {
      temp1 <- sum(pnorm( ga02%*%t(wij[1:n1,]) ,log.p=TRUE))
      #log(pnorm( ga02%*%t(wij[1:n1,]) )))
      temp2 <- sum(pnorm( ga02%*%t(wij[(n1+1):n,]) , lower.tail=FALSE, log.p=TRUE))
      #log(1-pnorm( ga02%*%t(wij[(n1+1):n,]) )))
      return(loglik <- -(temp1+temp2))
    }
  } else if(select.dist=="logistic") {
    fga <- function(ga02) {
      temp1 <- sum(plogis( ga02%*%t(wij[1:n1,]) ,log.p=TRUE) )
      #temp2 <- sum(log(1-plogis(ga02%*%t(wij[(n1+1):n,]))))
      #temp2 <- sum(plogis( -ga02%*%t(wij[(n1+1):n,]), log.p=TRUE ) )
      temp2 <- sum(plogis( ga02%*%t(wij[(n1+1):n,]) , lower.tail=FALSE, log.p=TRUE ) )
      return(loglik <- -(temp1+temp2))
    }
  }
  theta <- theta0

  ##ga_esti : (NRBin-ftn3.R)
  if (is.null(theta0)) # To give initial values
  {
    b0b1glm <- glm(yi[1:n1] ~ xij[1:n1,2:pdim], family=binomial(logit))
    ga_esti <- optim(c(rep(1,qdim)),fga)$par
    theta <- as.numeric(c(b0b1glm$coef,ga_esti)) # initial values)
  }

  ##ga_esti : (bin-funct.R)  (add observeY), don't use fga function
  #  if (is.null(theta0)) # To give initial values
  # {
  #    observeY <- !is.na(yi)
  #    b0b1glm <- glm(yi[1:n1] ~ xij[1:n1,2:pdim], family=binomial(logit))
  #    ga_esti <- glm(as.numeric(observeY) ~ wij[,2:qdim],family=binomial(logit)))
  #    theta <- as.numeric(c(b0b1glm$coef,ga_esti$coef)) # initial values)
  #  }


  alpha <- alpha0
  diff <- 1
  it <- 0; initmax <- itmax # To avoid an infinite loop
  likval <- like_bin(mi, xij, wij, theta, yi, alpha, select.dist)
  if(is.na(likval)) {
    return(list(theta<-NA, Jbar<-NA, likval<-NA)) # estimated parameters and Fisher observed information matrix
  }

  while ((diff > eps) && (it < itmax))
  {
    it <- it+1
    theta.old <- theta; likval.old <- likval

    # score and observed information matrix ------------------------------------------------------
    score_Jbar <- eval_score_obsinf_bin(mi, xij, wij, theta, yi, alpha, select.dist)

    score <- score_Jbar[[1]]
    Jbar <- score_Jbar[[2]]

    tCh <- tryCatch(is.positive.definite(Jbar),error=function(e){999})
    if(tCh != 999) {
      Jbar <- posdefify(Jbar) # To stabilize additionally (make it a pd matrix)
    } else {
      return(list(theta<-NA, Jbar<-NA, likval<-NA))
    }
    likval <- likval.old - 1; delta <- 1
    init <- 0
    while((likval < likval.old) && (init < initmax)) {

      theta <- theta.old + delta*solve(Jbar,score)
      init <- init+1

      #       ifelse(alpha < alpha0,likval <- -Inf,likval <- like(mi, xij, wij, theta, yi, alpha, negasso))

      likval <- like_bin(mi, xij, wij, theta, yi, alpha, select.dist)
      #       print(likval)
      if(is.na(likval)) likval <- -Inf

      delta <- delta/2
    }
    diff <- sum(abs(theta-theta.old))
  }
  list(theta<-theta, Jbar<-Jbar, likval<-likval) # estimated parameters and Fisher observed information matrix
}

#----------------------------------------------------------------------------------
like_bin <- function(mi, xij, wij, theta, yi, alpha, select.dist) {
  # log-likelihood function ---------------------------------------------
  pdim <- ncol(xij); qdim <- ncol(wij)
  pardim <- pdim+qdim
  n <- nrow(xij)
  n1 <- n-sum(is.nan(yi)) # count # of NaN's

  xibeta <- colSums(diag(theta[1:pdim])%*%t(xij))
  #xibeta[ which( xibeta <= -.Machine$double.xmax ) ] <- -.Machine$double.xmax
  #xibeta[ which( xibeta >= .Machine$double.xmax ) ] <- .Machine$double.xmax

  wigamma <- colSums(diag(theta[(pdim+1):(pdim+qdim)])%*%t(wij))

  #probi <- exp(xibeta)/(1+exp(xibeta))
  probi <- 1/(1+exp(-xibeta))
  probi[ which(probi < .Machine$double.eps) ] <- .Machine$double.eps
  probi[ which(probi > 1) ] <- 1 - .Machine$double.eps

  mui <- mi*probi

  alpha <- alpha

  taui <- wigamma
  taui[ which( taui <= -.Machine$double.xmax ) ] <- -.Machine$double.xmax
  taui[ which( taui >= .Machine$double.xmax ) ] <- .Machine$double.xmax

  etai <- alpha/mui
  etai[ which( etai <= -.Machine$double.xmax) ] <- -.Machine$double.xmax
  etai[ which( etai >= .Machine$double.xmax ) ] <- .Machine$double.xmax

  if(select.dist=="gumbel") {
    hyi <- exp(taui+etai*yi)
  } else if(select.dist=="normal") {
    hyi <- taui+etai*yi
    hyi[ which(hyi <= -.Machine$double.xmax) ] <- -.Machine$double.xmax
    hyi[ which(hyi >= .Machine$double.xmax) ] <- .Machine$double.xmax
  } else if(select.dist=="logistic") {
    hyi <- taui+etai*yi
    hyi[ which(hyi <= -.Machine$double.xmax) ] <- -.Machine$double.xmax
    hyi[ which(hyi >= .Machine$double.xmax) ] <- .Machine$double.xmax
  }

  #  if(!sum(hyi[1:n1]>0)==n1) {
  #    print("Give a larger lower bound for alpha0.")
  #    return(tempLik=NA)
  #  }

  if(select.dist=="gumbel"){
    #tempLik0 <- dbinom(yi[1:n1], mi[1:n1], probi[1:n1], log=TRUE) + pexp(hyi[1:n1], log.p=TRUE)
    tempLik0 <- dbinom(yi[1:n1], mi[1:n1], probi[1:n1]) * pexp(hyi[1:n1]) # add
    for ( i in 1:length(tempLik0)) {
      if (tempLik0[i] <= .Machine$double.eps) {
        tempLik0[i] <- .Machine$double.eps
      } else if (tempLik0[i] >= 1) {
        tempLik0[i] <- 1 - .Machine$double.eps
      }
    }
    tempLik0 <- log(tempLik0)
    #yi[1:n1]*log(probi[1:n1]/(1-probi[1:n1]))+ mi[1:n1]*log(1-probi[1:n1])+log(choose(mi[1:n1],yi[1:n1])) + log(1-exp(-hyi[1:n1]))
    tempLik1 <-  1- (1-probi[(n1+1):n])*pexp(exp(taui[(n1+1):n]))-probi[(n1+1):n]*pexp(exp(taui[(n1+1):n]+etai[(n1+1):n]))
    for ( i in 1:length(tempLik1)) {
      if (tempLik1[i] < 0) {
        tempLik1[i] <- .Machine$double.eps
      }
    }
    tempLik1 <- log(tempLik1)
  } else if(select.dist=="normal") {
    tempLik0 <- dbinom(yi[1:n1], mi[1:n1], probi[1:n1], log=TRUE) + pnorm(hyi[1:n1], log.p=TRUE)
    tempLik1 <- 1-(1-probi[(n1+1):n])*pnorm(taui[(n1+1):n])-probi[(n1+1):n]*pnorm(taui[(n1+1):n]+etai[(n1+1):n])
    for ( i in 1:length(tempLik1)) {
      if (tempLik1[i] < 0) {
        tempLik1[i] <- .Machine$double.eps
      }
    }
    tempLik1 <- log(tempLik1)
  } else if(select.dist=="logistic") {
    #tempLik0 <- dbinom(yi[1:n1], mi[1:n1], probi[1:n1], log=TRUE) + plogis(hyi[1:n1], log.p=TRUE)
    tempLik0 <- dbinom(yi[1:n1], mi[1:n1], probi[1:n1]) * plogis(hyi[1:n1]) # add
    for ( i in 1:length(tempLik0)) {
      if (tempLik0[i] < 0) {
        tempLik0[i] <- .Machine$double.eps
      }
    }
    tempLik0 <- log(tempLik0)

    tempLik1 <- 1-(1-probi[(n1+1):n])*plogis(taui[(n1+1):n]) - probi[(n1+1):n]*plogis(taui[(n1+1):n]+alpha/probi[(n1+1):n])
    for ( i in 1:length(tempLik1)) {
      if (tempLik1[i] < 0) {
        tempLik1[i] <- .Machine$double.eps
        }
    }
    tempLik1 <- log(tempLik1)
    #tempLik1 <- log(1-(1-probi[(n1+1):n])/(1+exp(-taui[(n1+1):n])) - probi[(n1+1):n]/(1+exp(-(taui[(n1+1):n] + alpha/(probi[(n1+1):n])))) )
  }
  tempLik <-  c(tempLik0, tempLik1)

  return(sum(tempLik))
}
