#' Fitting Sample Selection Models for Discrete Response Variables
#'
#' Function \code{DiSSMod} fits sample selection models for discrete random
#' variables, by suitably extending the formulation of the classical
#' Heckman model to the case of a discrete response, but retaining the
#' original conceptual framework. Maximum likelihood estimates are obtained
#' by Newton-Raphson iteration combined with use of profile likelihood.
#'
#' @details
#' The specification of the two linear models regulating the response variable and
#' the selection mechanism, as indicated in the \sQuote{Background} section,
#' is accomplished by two arguments of \code{formula} type,
#' denoted \code{response} and \code{selection}, respectively.
#' Each \code{formula} is specified with the same syntax of similar arguments in
#' standard functions such as \code{lm} and \code{glm}, with the restriction that
#' the intercept term (which is automatically included) must not be removed.
#'
#' The distributional assumptions associated to the \code{response} and \code{selection} components
#' are specified by the arguments \code{resp.dist} and \code{select.dist}, respectively.
#' Argument \code{select.dist} refers to the unobservable continuous variable of which we
#' observe only the dichotomous outcome Yes-No.
#'
#' In this respect, a remark is appropriate about the option \code{"Gumbel"} for \code{select.dist}.
#' This choice is equivalent to the adoption of an Exponential distribution of the selection variables
#' combined  an exponential transformation of the linear predictor of the
#' \code{selection} argument, as it is presented in Section 3.2 of Azzalini et al. (2019).
#' Also, it corresponds to work with the log-transformation of an Exponential variable,
#' which is essentially a Gumbel type of variable, up to a linear transformation with
#' respect to its more commonly employed parameterization.
#'
#' When \code{resp.dist} is \code{"Poisson"} or \code{"NegBinomial"} and \code{trunc.num} is missing,
#' a default choice is made; this equals \code{1.5*m} or \code{2*m} in the two respective cases,
#' where \code{m} denotes the maximum observed value of the response variable.
#'
#' Function \code{DiSSMOd} calls lower level functions, \code{nr.bin, nr.nbinom, nr.pois} and the others
#' for the actual numerical maximization of the log-likelihood via a Newton-Raphson iteration.
#'
#' Notice that the automatic initialization of the \code{alpha} search interval, when this argument is
#' missing, may change in future versions of the package.
#'
#' @concept generalized linear models
#' @concept selection models
#' @concept sample selection
#' @concept Heckman model
#' @concept discrete response
#' @concept maximum likelihood
#' @param response a formula for the response equation.
#' @param selection a formula for the selection equation.
#' @param data a data frame and data has to be included with the form of \code{data.frame}.
#' @param resp.dist a character for the distribution choice of the response variable,
#' \code{"bernoulli"} for Bernoulli distribution, \code{"poisson"} for Poisson distribution,
#' and \code{"negbinomial"} for Negative binomial distribution. Also, the character strings
#' can be abbreviated and can be upper or lower case as preferred.
#' @param select.dist a character for the distribution choice of the selection variable,
#' \code{"gumbel"} for Gumbel distribution, \code{"normal"} for Normal distribution,
#' and \code{"logistic"} for Logistic distribution. Also, the character strings
#' can be abbreviated and can be upper or lower case as preferred.
#' @param alpha a vector of \eqn{alpha} values on which the profile log-likelihood function is evaluated;
#' if the argument is missing, a set of values in the interval \code{(-10, 10)} is used for the initial search,
#' followed by a second search on a revised interval which depends on the outcome from the first search.
#' @param trunc.num an integer numeric constant used as the truncation point of an infine summation of probabilities
#' involved when \code{resp.dist} equals \code{"Poisson"} or \code{"NegBinomial"};
#' if the argument is missing, a default choice is made, as described in the \sQuote{Details} section. Notice: this
#' default choice of \code{trunc.num} may be subject to revision in some future version of the package,
#' and the argument \code{trunc.num} itselt may possibly be replaced by some other ingredient.
#' @param standard a logical value for the standardizing explanatory variables, if \code{TRUE} two types of values
#' (standardized and not) will be returned.
#' @param verbose an integer value for the level of printed details (values: 0|1|2); the default value is 1
#' which stands for shortly printed details. If the value is 2, more details are viewed such as
#' values of the log likelihood functions and iteration numbers. If the value is 0, there is no printed
#' detail.
#' @param eps a numeric value for the estimating parameters, which is needed for the step of the optimization.
#' If the sum of absolute differences between present step estimated parameters and former step
#' estimated parameters is smaller than \code{eps}, we assume that estimated parameters are
#' optimized.
#' @param itmax an integer stands for maximum number for the iteration of optimizing the parameters.
#' @return \code{DiSSMod} returns an object of class \code{"DiSSMod"},
#' which is a list containing following components:
#' @return
#' \item{call}{a matched call.}
#' \item{standard}{a logical value, stands for standardization or not.}
#' \item{st_loglik}{a vector containing the differences between log likelihoods and maximized log likelihood.}
#' \item{max_loglik}{a maximized log likelihood value.}
#' \item{mle_alpha}{a maximized likelihood estimator of alpha.}
#' \item{alpha}{a vector containing grids of the alpha}
#' \item{Nalpha}{a vector containing proper alpha, which does not have
#' \code{NA} value for corresponding log likelihood.}
#' \item{num_NA}{a number of \code{NA} values of log likelihoods.}
#' \item{n_select}{a number of selected response variables.}
#' \item{n_all}{a number of all response variables.}
#' \item{estimate_response}{estimated values for the response model.}
#' \item{std_error_response}{estimated standard errors for the response model.}
#' \item{estimate_selection}{estimated values for the selection model.}
#' \item{std_error_selection}{estimated standard errors for the selection model.}
#' @importFrom MASS glm.nb
#' @importFrom sfsmisc posdefify
#' @importFrom matrixcalc is.positive.definite
#' @importFrom grDevices extendrange
#' @importFrom graphics abline plot points rug
#' @importFrom stats approxfun binomial contrasts dbinom
#' dexp dlogis dnbinom dnorm dpois glm model.matrix
#' model.response na.pass optim pexp plogis pnorm
#' poisson qchisq qnorm sd symnum uniroot
#' @import psych
#' @examples
#' set.seed(45)
#' data(DoctorRWM, package = "DiSSMod")
#' n0 <- 600
#' set.n0 <- sample(1:nrow(DoctorRWM), n0)
#' reduce_DoctorRWM <- DoctorRWM[set.n0,]
#' result0 <- DiSSMod(response = as.numeric(DOCVIS > 0) ~ AGE + INCOME_SCALE + HHKIDS + EDUC + MARRIED,
#'                    selection = PUBLIC ~ AGE + EDUC + FEMALE,
#'                    data = reduce_DoctorRWM, resp.dist="bernoulli", select.dist = "normal",
#'                    alpha = seq(-5.5, -0.5, length.out = 21), standard = TRUE)
#'
#' print(result0)
#'
#' data(CreditMDR, package = "DiSSMod")
#' n1 <- 600
#' set.n1 <- sample(1:nrow(CreditMDR), n1)
#' reduce_CreditMDR <- CreditMDR[set.n1,]
#' result1 <- DiSSMod(response = MAJORDRG ~ AGE + INCOME + EXP_INC,
#'                    selection = CARDHLDR ~ AGE + INCOME + OWNRENT + ADEPCNT + SELFEMPL,
#'                    data = reduce_CreditMDR, resp.dist="poi", select.dist = "logis",
#'                    alpha = seq(-0.3, 0.3,length.out = 21), standard = FALSE, verbose = 1)
#'
#' print(result1)
#'
#' @references Azzalini, A., Kim, H.-M. and Kim, H.-J. (2019) Sample selection
#' models for discrete and other non-Gaussian response variables.
#'  \emph{Statistical Methods & Applications}, \strong{28}, 27--56. First online 30 March 2018.
#' \url{https://doi.org/10.1007/s10260-018-0427-1}
#' @seealso The functions \code{\link[DiSSMod]{summary.DiSSMod}}, \code{\link[DiSSMod]{coef.DiSSMod}},
#' \code{\link[DiSSMod]{confint.DiSSMod}}, \code{\link[DiSSMod]{plot.DiSSMod}}
#'  are used to obtain and print a summary, coefficients, confidence interval and
#'  plot of the results.
#' @seealso The generic function \code{\link[stats]{logLik}} is used to obtain maximum log likelihood of the
#' result.
#' @seealso See also \code{\link[stats]{lm}}, \code{\link[stats]{glm}} and
#' \code{\link[stats]{formula}}.
#' @section Background:
#' Function \code{DiSSMod} fits sample selection models for discrete random variables,
#' by suitably extending the formulation of the classical Heckman model to the case of a discrete response,
#' but retaining the original conceptual framework.
#' This logic involves the following key ingredients: (1) a linear model indicating which explanatory variables
#' influence the response variable; (2) a linear model indicating  which (possibly different) explanatory variables,
#' besides the response variable itself, influence  a `selection variable', which is intrinsically continuous but
#' we only observe a dichotomous outcome from it, of type Yes-No, which selects which are the observed response cases;
#' (3) distributional assumptions on the response and the selection variable.
#'
#' The data fitting method is maximum likelihood estimation (MLE), which operates in two steps:
#' (i) for each given value of parameter \eqn{alpha} which regulates the level of selection,
#'  MLE is performed for all the remaining parameters, using a Newton-Raphson iteration;
#' (ii) a scan of the \eqn{alpha} axis builds the  profile log-likelihood function and
#'  its maximum point represents the overall MLE.
#'
#' A detailed account of the underlying theory and the operational methodology is provided by Azzalini et al. (2019).
#'
#' @export
DiSSMod <- function(response, selection, data, resp.dist, select.dist, alpha,
                    trunc.num, standard = FALSE, verbose = 1, eps = 1e-7,
                    itmax = 1e+3)
{
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)

  m1 <- match(c("response", "data", "subset", "offset"), names(mf), 0)
  mf1 <- mf[c(1, m1)]
  mf1$drop.unused.levels <- TRUE
  mf1[[1]] <- quote(model.frame)
  names(mf1)[2] <- "formula"
  mf1 <- eval(mf1, parent.frame())
  mt1 <- attr(mf1, "terms")

  m2 <- match(c("selection", "data", "subset"), names(mf), 0)
  mf2 <- mf[c(1, m2)]
  mf2$drop.unused.levels <- TRUE
  mf2$na.action <- na.pass
  names(mf2)[2] <- "formula"
  mf2[[1]] <- quote(model.frame)
  mf2 <- eval(mf2, parent.frame())
  mt2 <- attr(mf2, "terms")

  resp.dist <- match.arg(tolower(resp.dist), c("bernoulli", "poisson", "negbinomial"))
  select.dist <- match.arg(tolower(select.dist), c("gumbel", "normal", "logistic"))

  pre.search <- missing(alpha)

  y.out <- model.response(mf1, "numeric")
  u.sel <- model.response(mf2, "numeric")

  m <- max(y.out, na.rm = TRUE)
  if(missing(trunc.num)) {
    if(resp.dist=="poisson") {
      trunc.num <- 1.5*m
    } else if(resp.dist=="negbinomial") {
      trunc.num <- 2*m
    } else {
      trunc.num <- NULL }
  }

  n <- length(y.out)
  n1 <- sum(u.sel==1)

  yi <- y.out[u.sel==1]
  yi <- c(yi,rep(NaN,n-n1))

  x.mat <- model.matrix(mt1, mf1, contrasts)[,-1]
  w.mat <- model.matrix(mt2, mf2, contrasts)[,-1]

  x.name <- names(as.data.frame(x.mat))
  w.name <- names(as.data.frame(w.mat))

  x.mat <- apply(as.matrix(x.mat), 2, function(x, a, b){c(x[a],x[b])}, a = c(u.sel==1), b = c(u.sel==0) )
  w.mat <- apply(as.matrix(w.mat), 2, function(x, a, b){c(x[a],x[b])}, a = c(u.sel==1), b = c(u.sel==0) )

  if(standard){
    mean.x.vec <- apply(x.mat, 2, mean)
    mean.w.vec <- apply(w.mat, 2, mean)
    sd.x.vec <- apply(x.mat, 2, sd)
    sd.w.vec <- apply(w.mat, 2, sd)

    x.mat <- apply(as.matrix(x.mat), 2, function(x){(x-mean(x))/sd(x)} )
    w.mat <- apply(as.matrix(w.mat), 2, function(x){(x-mean(x))/sd(x)} )

    cm <- diag(1/c(1, sd.x.vec, 1, sd.w.vec))
    cm[1, c(2:(ncol(x.mat)+1))] <- c(-mean.x.vec/sd.x.vec)
    cm[(ncol(x.mat)+2), c((ncol(x.mat)+3):(ncol(x.mat)+ncol(w.mat)+2))] <- c(-mean.w.vec/sd.w.vec)
  }

  one <- rep(1,n)

  xij <- cbind(one,x.mat)
  wij <- cbind(one,w.mat)

  if(pre.search) {
    alpha <- seq(-10, 10, length=21)
    if(verbose > 0) cat("Search preliminary alpha interval (", format(range(alpha)), ")\n")
  }

  leng <- length(alpha)

  loglik <- rep(0,leng)
  theta_Finf <- NULL
  mi <- rep(1,n)

  if(verbose == 1) cat("Running up to", leng, ": ")
  for (i in 1:leng) {
    alpha0 <- alpha[i]
    theta0 <- NULL
    temp <- switch(resp.dist,
                   "bernoulli" = nr.bin(mi, xij, wij, theta0, yi, alpha0, select.dist, eps, itmax),
                   "poisson" = nr.pois(xij, wij, theta0, yi, alpha0, select.dist, eps, itmax, trunc.num),
                   "negbinomial" = nr.nbinom(xij, wij, theta0, yi, alpha0, select.dist, eps, itmax, trunc.num)
    )
    loglik[i] <- temp[[3]]
    theta_Finf <- c(theta_Finf,temp)
    if(verbose == 1) cat(i, "")
    if(verbose > 1) cat("i, alpha[i], logLik[i]:", i, format(alpha0), format(loglik[i]), "\n")
  }
  if(verbose == 1) cat("\n")
  numind <- which.max(loglik)
  if(pre.search) {
    alpha. <- alpha[numind]
    bad <- is.na(loglik) | (loglik[numind] - loglik > 15)
    shift <- loglik[numind] - qchisq(0.9, 1)
    if(sum(!bad) < 3) {
      message("This log-likelihood function takes on peculiar values.")
      if(verbose < 2) print(rbind(alpha, logLik=loglik))
      stop("I give up. Sorry.")
    }
    af <- approxfun(alpha[!bad], loglik[!bad] - shift, method="linear")
    vaf <- Vectorize(af)
    alpha1 <- try(uniroot(vaf, c(min(alpha[!bad]), alpha.)), silent = TRUE)
    alpha1 <- if(class(alpha1)=="list") alpha1$root else min(alpha)
    alpha2 <- try(uniroot(vaf, c(alpha., max(alpha[!bad]))), silent = TRUE)
    alpha2 <- if(class(alpha2)=="list") alpha2$root else max(alpha)
    r <- extendrange(c(alpha1, alpha2), f=0.2)
    leng <- 26
    alpha <- seq(r[1], r[2], length=leng)
    loglik <- rep(0,leng)
    theta_Finf <- NULL
    mi <- rep(1,n)
    if(verbose > 0) cat("Search in revised alpha interval (", format(range(alpha)), ")\n")
    if(verbose == 1) cat("Running up to", leng, ": ")
    for (i in 1:leng) {
      alpha0 <- alpha[i]
      theta0 <- NULL
      temp <- switch(resp.dist,
                     "bernoulli" = nr.bin(mi, xij, wij, theta0, yi, alpha0, select.dist, eps, itmax),
                     "poisson" = nr.pois(xij, wij, theta0, yi, alpha0, select.dist, eps, itmax, trunc.num),
                     "negbinomial" = nr.nbinom(xij, wij, theta0, yi, alpha0, select.dist, eps, itmax, trunc.num)
      )
      loglik[i] <- temp[[3]]
      theta_Finf <- c(theta_Finf,temp)
      if(verbose == 1) cat(i, "")
      if(verbose > 1) cat("i, alpha[i], logLik[i]:", i, format(alpha0), format(loglik[i]), "\n")
    }
    if(verbose == 1) cat("\n")
    numind <- which.max(loglik)
  }

  estPar <- theta_Finf[[(numind-1)*3+1]]
  FisherInf <- theta_Finf[[(numind-1)*3+2]]
  if(resp.dist=="negbinomial"){
    FisherInf_b <- FisherInf
    stderr_b <- sqrt(diag(solve(  FisherInf_b  )))[ncol(x.mat)+ncol(w.mat)+3]
    psi <- estPar[ncol(x.mat)+ncol(w.mat)+3] <- exp(estPar[ncol(x.mat)+ncol(w.mat)+3])
    inv.b.mat <- diag(c(rep(1, ncol(x.mat)+ncol(w.mat)+2), 1/psi))
    FisherInf <- inv.b.mat %*% FisherInf %*% inv.b.mat
  }
  stError <- sqrt(diag(solve(  FisherInf  )))

  if(standard) {
    ori_estPar <- as.vector(cm %*% estPar)
    cov_mle <- solve(FisherInf)
    cov_MLE <- cm %*% cov_mle %*% t(cm)
    ori_stError <- sqrt(diag(cov_MLE))

    if(resp.dist == "negbinomial"){
      cov_mle_b <- solve(FisherInf_b)
      cov_MLE_b <- cm %*% cov_mle_b %*% t(cm)
      ori_stderr_b <- sqrt(diag(cov_MLE_b))[ncol(x.mat)+ncol(w.mat)+3]
    }

    ori_estimate_response <- ori_estPar[1:(ncol(x.mat)+1)]
    ori_std_error_response <- ori_stError[1:(ncol(x.mat)+1)]
    ori_estimate_selection <- ori_estPar[(ncol(x.mat)+2):(ncol(x.mat)+ncol(w.mat)+2)]
    ori_std_error_selection <- ori_stError[(ncol(x.mat)+2):(ncol(x.mat)+ncol(w.mat)+2)]
  }

  if(resp.dist=="negbinomial") {
    theta <- estPar[ncol(x.mat)+ncol(w.mat)+3]
    theta.stderr <- stError[ncol(x.mat)+ncol(w.mat)+3]
    if(standard) {
      ori_theta <- ori_estPar[ncol(x.mat)+ncol(w.mat)+3]
      ori_theta.stderr <- ori_stError[ncol(x.mat)+ncol(w.mat)+3]
    }
  }

  Nloglik <- loglik[!is.na(loglik)]
  Nalpha <- alpha[!is.na(loglik)]
  numNA <- sum(is.na(loglik))

  maxloglik <- max(Nloglik)
  mle_alpha <- Nalpha[which.max(Nloglik)]

  stloglik <- loglik - maxloglik

  result <- list( call = cl, standard = standard, st_loglik = stloglik,
                  max_loglik = maxloglik, mle_alpha = mle_alpha,
                  alpha = alpha, Nalpha = Nalpha, num_NA = numNA,
                  n_select = n1, n_all = n,
                  estimate_response = estPar[1:(ncol(x.mat)+1)],
                  std_error_response = stError[1:(ncol(x.mat)+1)],
                  estimate_selection = estPar[(ncol(x.mat)+2):(ncol(x.mat)+ncol(w.mat)+2)],
                  std_error_selection = stError[(ncol(x.mat)+2):(ncol(x.mat)+ncol(w.mat)+2)]
  )
  if(resp.dist=="negbinomial"){
    result$theta <- theta
    result$theta_std_error <- theta.stderr
    result$std_error_b <- stderr_b
  }
  if(standard){
    result$ori_estimate_response <- ori_estimate_response
    result$ori_std_error_response <- ori_std_error_response
    result$ori_estimate_selection <- ori_estimate_selection
    result$ori_std_error_selection <- ori_std_error_selection
    if(resp.dist=="negbinomial"){
      result$ori_theta <- ori_theta
      result$ori_theta_std_error <- ori_theta.stderr
      result$ori_b_std_error <- ori_stderr_b
    }
  }
  result$names_response <- c( "(Intercept)", x.name )
  result$names_selection <- c( "(Intercept)", w.name )

  class(result) <- "DiSSMod"
  result
}

#' @method print DiSSMod
#' @export
print.DiSSMod <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\nCall:\n")
  print(x$call)
  cat("\nMLE of alpha:", x$mle_alpha, "\n")
  cat("Maximum log likelihood:", x$max_loglik, "\n")
  if(x$standard) cat("\n--------------- Standardized Result ---------------\n")
  cat("\nResponse coefficients:\n")
  print(t(data.frame(Estimate = x$estimate_response,
                     row.names = x$names_response, check.names=FALSE)), digits = digits )
  cat("\nSelection coefficients:\n")
  print(t(data.frame(Estimate = x$estimate_selection,
                     row.names = x$names_selection, check.names=FALSE)), digits = digits )
  cat("\n")
  if(x$standard) {
    cat("--------------- Original Scale Result ---------------\n")
    cat("\nResponse coefficients:\n")
    print(t(data.frame(Estimate = x$ori_estimate_response,
                       row.names = x$names_response, check.names=FALSE)), digits = digits )
    cat("\nSelection coefficients:\n")
    print(t(data.frame(Estimate = x$ori_estimate_selection,
                       row.names = x$names_selection, check.names=FALSE)), digits = digits )
    cat("\n")
  }
}

#' Summarizing Discrete Sample Selection Model Fits
#'
#' \code{summary} method for a class \code{"DiSSMod"}.
#'
#' @param object an object of class \code{"DiSSMod"} made by the function \code{DiSSMod}.
#' @param x an object of class \code{"summary.DiSSMod"}.
#' @param digits a numeric number of significant digits.
#' @param \dots additional control argument is as follows.
#' \itemize{
#' \item \code{level}: an option for controlling the significance level of confidence interval.
#' It has to be given in probability between 0 and 1. Initial level is set to \eqn{1 - \alpha = 0.95}.
#' }
#' @details If \code{standard} equals \code{TRUE}, \code{summary} also additionally returns summary
#' statistics of standardized results. Otherwise, it just returns summary statistics as similar statistics
#' as the generic function \code{summary}.
#' @return The function \code{summary.DiSSMod} returns a list of summary statistics of the fitted
#' discrete sample selection model given in \code{object}.
#' @return The components, which are not duplicated from the \code{object}, are as follows:
#' @return
#' \item{z.value_response}{Z statistics (normal distribution) for coefficients of response equation.}
#' \item{z.value_selection}{Z statistics (normal distribution) for coefficients of selection equation.}
#' \item{CI_alpha}{confidence interval of the parameter \code{alpha}.}
#' \item{level}{a numeric value between 0 and 1 for controlling the significance level of confidence interval.
#'  Initial level is set to \eqn{1 - \alpha = 0.95}.}
#' @examples
#' # example continued from DiSSMod
#' set.seed(45)
#' data(DoctorRWM, package = "DiSSMod")
#' n0 <- 600
#' set.n0 <- sample(1:nrow(DoctorRWM), n0)
#' reduce_DoctorRWM <- DoctorRWM[set.n0,]
#' result0 <- DiSSMod(response = as.numeric(DOCVIS > 0) ~ AGE + INCOME_SCALE + HHKIDS + EDUC + MARRIED,
#'                    selection = PUBLIC ~ AGE + EDUC + FEMALE,
#'                    data = reduce_DoctorRWM, resp.dist="bernoulli", select.dist = "normal",
#'                    alpha = seq(-5.5, -0.5, length.out = 21), standard = TRUE)
#'
#' summary(result0, level = 0.90)
#'
#' data(CreditMDR, package = "DiSSMod")
#' n1 <- 600
#' set.n1 <- sample(1:nrow(CreditMDR), n1)
#' reduce_CreditMDR <- CreditMDR[set.n1,]
#' result1 <- DiSSMod(response = MAJORDRG ~ AGE + INCOME + EXP_INC,
#'                    selection = CARDHLDR ~ AGE + INCOME + OWNRENT + ADEPCNT + SELFEMPL,
#'                    data = reduce_CreditMDR, resp.dist="poi", select.dist = "logis",
#'                    alpha = seq(-0.3, 0.3,length.out = 21), standard = FALSE, verbose = 1)
#'
#' summary(result1)
#'
#' @seealso See also \code{\link{DiSSMod}} and \code{\link[base]{summary}}.
#' @method summary DiSSMod
#' @export
summary.DiSSMod <- function(object, ...)
{
  z.value_response <- object$estimate_response/object$std_error_response
  z.value_selection <- object$estimate_selection/object$std_error_selection
  if(object$standard){
    ori_z.value_response <- object$ori_estimate_response/object$ori_std_error_response
    ori_z.value_selection <- object$ori_estimate_selection/object$ori_std_error_selection
  }
  ci <- confint.DiSSMod(object, ...)
  result <- list( call = object$call, standard = object$standard, n_all = object$n_all, n_select = object$n_select,
                  max_loglik = object$max_loglik, mle_alpha = object$mle_alpha,
                  estimate_response = object$estimate_response, std_error_response = object$std_error_response,
                  estimate_selection = object$estimate_selection, std_error_selection = object$std_error_selection,
                  z.value_response = z.value_response, z.value_selection = z.value_selection
  )
  if(object$standard){
    result$ori_estimate_response <- object$ori_estimate_response
    result$ori_std_error_response <- object$ori_std_error_response
    result$ori_estimate_selection <- object$ori_estimate_selection
    result$ori_std_error_selection <- object$ori_std_error_selection
    result$ori_z.value_response <- ori_z.value_response
    result$ori_z.value_selection <- ori_z.value_selection
  }
  if(!is.null(object$theta)){
    result$theta <- object$theta
    result$theta_std_error <- object$theta_std_error
    if(object$standard){
      result$ori_theta <- object$ori_theta
      result$ori_theta_std_error <- object$ori_theta_std_error
    }
  }
  result$CI_alpha <- ci$alpha
  result$level <- ci$level
  result$names_response <- object$names_response
  result$names_selection <- object$names_selection

  class(result) <- "summary.DiSSMod"
  result
}

#' @rdname summary.DiSSMod
#' @method print summary.DiSSMod
#' @export
print.summary.DiSSMod <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\nCall:\n")
  print(x$call)
  cat("\nNumber of observations:", x$n_all ,", selected (1):",
      x$n_select,", not selected (0):", x$n_all - x$n_select,"\n")
  cat("Maximum log likelihood:", format(x$max_loglik, nsmall = 2))
  cat("\nMLE of alpha:",x$mle_alpha,"\n")
  cat("\n", x$level*100, "% Confidence interval of alpha:\n")
  cat("\n", format(x$CI_alpha, digits=digits), "\n")
  if(x$standard) cat("\n--------------- Standardized Result ---------------\n")
  cat("\nResponse equation:\n")
  out <- data.frame( estimate = x$estimate_response, std.error = x$std_error_response,
                     z.value = x$z.value_response, pval = 2*pnorm(abs(x$z.value_response), lower.tail = F) )
  Signif.out <- symnum(out$pval, corr = FALSE, na = FALSE,
                       cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                       symbols = c("***", "**", "*", ".", " "))
  out <- cbind(out, format(Signif.out))
  row.names(out) <- x$names_response ; colnames(out) <- c("Estimate", "Std. Error", "Z Value", "Pr(>|Z|)", "")
  print(out, digits = digits)
  if(!is.null(x$theta)){
    cat("\n              Theta: ", x$theta)
    cat("\n          Std. Err.: ", x$theta_std_error, "\n")
  }
  cat("\nSelection equation:\n")
  sel <- data.frame( estimate = x$estimate_selection, std.error = x$std_error_selection,
                     z.value = x$z.value_selection, pval = 2*pnorm(abs(x$z.value_selection), lower.tail = F) )
  Signif.sel <- symnum(sel$pval, corr = FALSE, na = FALSE,
                       cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                       symbols = c("***", "**", "*", ".", " "))
  sel <- cbind(sel, format(Signif.sel))
  row.names(sel) <- x$names_selection ; colnames(sel) <- c("Estimate", "Std. Error", "Z Value", "Pr(>|Z|)", "")
  print(sel, digits = digits)

  cat("---\n")
  cat("Signif. codes:", 0, "'***'", 0.001, "'**'", 0.01, "'*'", 0.05, "'.'", 0.1, "' '", 1, "\n\n")
  if(x$standard) {
    cat("--------------- Original Scale Result ---------------\n")
    cat("\nResponse equation:\n")
    ori_out <- data.frame( estimate = x$ori_estimate_response, std.error = x$ori_std_error_response,
                           z.value = x$ori_z.value_response, pval = 2*pnorm(abs(x$ori_z.value_response), lower.tail = F) )
    ori_Signif.out <- symnum(ori_out$pval, corr = FALSE, na = FALSE,
                             cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                             symbols = c("***", "**", "*", ".", " "))
    ori_out <- cbind(ori_out, format(ori_Signif.out))
    row.names(ori_out) <- x$names_response ; colnames(ori_out) <- c("Estimate", "Std. Error", "Z Value", "Pr(>|Z|)", "")
    print(ori_out, digits = digits)
    if(!is.null(x$ori_theta)){
      cat("\n              Theta: ", x$ori_theta)
      cat("\n          Std. Err.: ", x$ori_theta_std_error, "\n")
    }
    cat("\nSelection equation:\n")
    ori_sel <- data.frame( estimate = x$ori_estimate_selection, std.error = x$ori_std_error_selection,
                           z.value = x$ori_z.value_selection, pval = 2*pnorm(abs(x$ori_z.value_selection), lower.tail = F) )
    ori_Signif.sel <- symnum(ori_sel$pval, corr = FALSE, na = FALSE,
                             cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                             symbols = c("***", "**", "*", ".", " "))
    ori_sel <- cbind(ori_sel, ori_Signif.sel)
    row.names(ori_sel) <- x$names_selection ; colnames(ori_sel) <- c("Estimate", "Std. Error", "Z Value", "Pr(>|Z|)", "")
    print(ori_sel, digits = digits)

    cat("---\n")
    cat("Signif. codes:", 0, "'***'", 0.001, "'**'", 0.01, "'*'", 0.05, "'.'", 0.1, "' '", 1, "\n\n")
  }
}

#' @method logLik DiSSMod
#' @export
logLik.DiSSMod <- function(object, ...)
{
  result <- object$max_loglik
  attr(result, "df") <- 1 + length(c(object$estimate_response, object$estimate_response, object$theta))
  class(result) <- "logLik"
  result
}

#' Getting Coefficients of Discrete Sample Selection Model Fits
#'
#' \code{coef} method for a class \code{"DiSSMod"}.
#'
#' @param object an object of class "DiSSMod" made by the function \code{DiSSMod}.
#' @param only a character value for choosing specific variable's coefficients. Initial value is
#' @param \dots not used, but exists because of the compatibility.
#' \code{NULL}, which shows all variable's coefficients. If \code{"response"} is written,
#'  only coefficients for response variables will be returned, and if \code{"selection"} is written,
#'  only coefficients for selection variables will be returned.
#' @return a numeric vector or a list is given.
#' @details It looks as similar as the generic function \code{coef}, but this case there
#' are two equations. Therefore, there exist little differences.
#' @examples
#' # example continued from DiSSMod
#' set.seed(45)
#' data(DoctorRWM, package = "DiSSMod")
#' n0 <- 600
#' set.n0 <- sample(1:nrow(DoctorRWM), n0)
#' reduce_DoctorRWM <- DoctorRWM[set.n0,]
#' result0 <- DiSSMod(response = as.numeric(DOCVIS > 0) ~ AGE + INCOME_SCALE + HHKIDS + EDUC + MARRIED,
#'                    selection = PUBLIC ~ AGE + EDUC + FEMALE,
#'                    data = reduce_DoctorRWM, resp.dist="bernoulli", select.dist = "normal",
#'                    alpha = seq(-5.5, -0.5, length.out = 21), standard = TRUE)
#'
#' coef(result0)
#'
#' data(CreditMDR, package = "DiSSMod")
#' n1 <- 600
#' set.n1 <- sample(1:nrow(CreditMDR), n1)
#' reduce_CreditMDR <- CreditMDR[set.n1,]
#' result1 <- DiSSMod(response = MAJORDRG ~ AGE + INCOME + EXP_INC,
#'                    selection = CARDHLDR ~ AGE + INCOME + OWNRENT + ADEPCNT + SELFEMPL,
#'                    data = reduce_CreditMDR, resp.dist="poi", select.dist = "logis",
#'                    alpha = seq(-0.3, 0.3,length.out = 21), standard = FALSE, verbose = 1)
#'
#' coef(result1)
#'
#' @seealso See also \code{\link{DiSSMod}} and \code{\link[stats]{coef}}.
#' @method coef DiSSMod
#' @export
coef.DiSSMod <- function(object, only = NULL, ...)
{
  if(object$standard){
    if (is.null(only)){
      names_response <- object$names_response
      names_selection <- object$names_selection
      result <- list( response = object$estimate_response,
                      selection = object$estimate_selection,
                      original_response = object$ori_estimate_response,
                      original_selection = object$ori_estimate_selection
      )
      names(result$response) <- names_response
      names(result$original_response) <- names_response
      names(result$selection) <- names_selection
      names(result$original_selection) <- names_selection
    } else if (only=="response"){
      names_response <- object$names_response
      result <- list( response = object$estimate_response,
                      original_response = object$ori_estimate_response
      )
      names(result$response) <- names_response
      names(result$original_response) <- names_response
    } else if (only=="selection"){
      names_selection <- object$names_selection
      result <- list( selection <- object$estimate_selection,
                      original_selection <- object$ori_estimate_selection
      )
      names(result$selection) <- names_selection
      names(result$original_selection) <- names_selection
    } else {
      stop("'only' has to be defined properly")
    }
  } else {
    if (is.null(only)){
      names_response <- object$names_response
      names_selection <- object$names_selection
      result <- list( response = object$estimate_response,
                      selection = object$estimate_selection
      )
      names(result$response) <- names_response
      names(result$selection) <- names_selection
    } else if (only=="response"){
      names_response <- object$names_response
      result <- object$estimate_response
      names(result) <- names_response
    } else if (only=="selection"){
      names_selection <- object$names_selection
      result <- object$estimate_selection
      names(result) <- names_selection
    } else {
      stop("'only' has to be defined properly")
    }
  }
  result
}

#' Getting Confidence Intervals for Parameters of Discrete Sample Selection Model Fits
#'
#' \code{confint} method for a class \code{"DiSSMod"}.
#'
#' @param object an object of class "DiSSMod" made by the function \code{DiSSMod}.
#' @param parm not used, but it exists for compatibility reasons.
#' @param level a numeric value between 0 and 1 for controlling the significance level of confidence interval;
#' default value is 0.95.
#' @param \dots not used, but it exists for compatibility reasons.
#' @return a list, containing \code{level} and \code{confidence intervals} for parameters, is given.
#'
#' @method confint DiSSMod
#' @examples
#' # example continued from DiSSMod
#' set.seed(45)
#' data(DoctorRWM, package = "DiSSMod")
#' n0 <- 600
#' set.n0 <- sample(1:nrow(DoctorRWM), n0)
#' reduce_DoctorRWM <- DoctorRWM[set.n0,]
#' result0 <- DiSSMod(response = as.numeric(DOCVIS > 0) ~ AGE + INCOME_SCALE + HHKIDS + EDUC + MARRIED,
#'                    selection = PUBLIC ~ AGE + EDUC + FEMALE,
#'                    data = reduce_DoctorRWM, resp.dist="bernoulli", select.dist = "normal",
#'                    alpha = seq(-5.5, -0.5, length.out = 21), standard = TRUE)
#'
#' confint(result0, level = 0.90)
#'
#' data(CreditMDR, package = "DiSSMod")
#' n1 <- 600
#' set.n1 <- sample(1:nrow(CreditMDR), n1)
#' reduce_CreditMDR <- CreditMDR[set.n1,]
#' result1 <- DiSSMod(response = MAJORDRG ~ AGE + INCOME + EXP_INC,
#'                    selection = CARDHLDR ~ AGE + INCOME + OWNRENT + ADEPCNT + SELFEMPL,
#'                    data = reduce_CreditMDR, resp.dist="poi", select.dist = "logis",
#'                    alpha = seq(-0.3, 0.3,length.out = 21), standard = FALSE, verbose = 1)
#'
#' confint(result1)
#'
#' @seealso See also \code{\link[stats]{confint}}, \code{\link[DiSSMod]{DiSSMod}} and \code{\link[DiSSMod]{summary.DiSSMod}}.
#' @export
confint.DiSSMod <- function(object, parm, level = 0.95, ...)
{
  coef <- coef(object)
  cf_resp <- coef$response
  cf_select <- coef$selection
  if(object$standard){
    cf_ori_resp <- coef$original_response
    cf_ori_select <- coef$original_selection
  }
  pnames_resp <- names(cf_resp)
  pnames_select <- names(cf_select)
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  fac <- qnorm(a)
  pct <- paste(as.character(a*100), "%")
  ci_resp <- array( NA, dim = c(length(pnames_resp), 2), dimnames = list(pnames_resp, pct) )
  ci_select <- array( NA, dim = c(length(pnames_select), 2), dimnames = list(pnames_select, pct) )
  if(object$standard){
    ci_ori_resp <- array( NA, dim = c(length(pnames_resp), 2), dimnames = list(pnames_resp, pct) )
    ci_ori_select <- array( NA, dim = c(length(pnames_select), 2), dimnames = list(pnames_select, pct) )
  }
  ses_resp <- object$std_error_response
  ses_select <- object$std_error_selection
  if(object$standard){
    ses_ori_resp <- object$ori_std_error_response
    ses_ori_select <- object$ori_std_error_selection
  }
  ci_resp[] <- cf_resp[pnames_resp] + ses_resp %o% fac
  ci_select[] <- cf_select[pnames_select] + ses_select %o% fac
  if(object$standard){
    ci_ori_resp[] <- cf_ori_resp[pnames_resp] + ses_ori_resp %o% fac
    ci_ori_select[] <- cf_ori_select[pnames_select] + ses_ori_select %o% fac
  }
  q <- qchisq(level, 1)
  na <- is.na(object$st_loglik)
  af <- approxfun(object$Nalpha[!na], 2*object$st_loglik[!na] + q)
  vaf <- Vectorize(af)
  alpha1 <- try(uniroot(vaf, c(min(object$Nalpha), object$mle_alpha)), silent = TRUE)
  alpha1 <- if(class(alpha1)=="list") alpha1$root else NA
  alpha2 <- try(uniroot(vaf, c(object$mle_alpha, max(object$Nalpha))), silent = TRUE)
  alpha2 <- if(class(alpha2)=="list") alpha2$root else NA
  alpha <- c(alpha1, alpha2)
  names(alpha) <- pct
  if(!is.null(object$theta)){
    theta <- c(log(object$theta) + (object$std_error_b * fac))
    theta <- exp(theta)
    names(theta) <- pct
    if(object$standard){
      ori_theta <- c(log(object$ori_theta) + (object$ori_std_error_b * fac))
      ori_theta <- exp(ori_theta)
      names(ori_theta) <- pct
    }
  }

  if(object$standard){
    result <- list(level = level, alpha = alpha, response = ci_resp, selection = ci_select,
                   original_response = ci_ori_resp, original_selection = ci_ori_select)
    if(!is.null(object$theta)){
      result$theta <- theta
      result$ori_theta <- ori_theta
    }
  } else {
    result <- list(level = level, alpha = alpha, response = ci_resp, selection = ci_select)
    if(!is.null(object$theta)){
      result$theta <- theta
    }
  }
  result
}

#' Relative Log Likelihood Plot for Discrete Sample Selection Model Fits
#'
#' \code{plot} method for a class \code{"DiSSMod"}.
#'
#' @method plot DiSSMod
#' @param x an object of class "DiSSMod" made by the function \code{DiSSMod}.
#' @param \dots additional control argument is as follows.
#' \itemize{
#' \item \code{level}: an option for controlling the significance level of confidence interval.
#' It has to be given in probability between 0 and 1. Initial level is set to \eqn{1 - \alpha = 0.95}.
#' }
#' @details Function \code{plot} draws a convex line due to the values of twice relative
#' log likelihoods by using the profile likelihood approach with
#' following the grids of \code{alpha}. If confidence interval created from
#' the function \code{confint} exists between the maximum and minimum value of the \code{alpha},
#' there will be two points drawn with the color red. Also, the Maximum Likelihood Estimator (MLE)
#' of \code{alpha} can be seen easily, if it exists between the maximum and minimum value of
#' the \code{alpha}.
#' @examples
#' # example continued from DiSSMod
#' set.seed(45)
#' data(DoctorRWM, package = "DiSSMod")
#' n0 <- 600
#' set.n0 <- sample(1:nrow(DoctorRWM), n0)
#' reduce_DoctorRWM <- DoctorRWM[set.n0,]
#' result0 <- DiSSMod(response = as.numeric(DOCVIS > 0) ~ AGE + INCOME_SCALE + HHKIDS + EDUC + MARRIED,
#'                    selection = PUBLIC ~ AGE + EDUC + FEMALE,
#'                    data = reduce_DoctorRWM, resp.dist="bernoulli", select.dist = "normal",
#'                    alpha = seq(-5.5, -0.5, length.out = 21), standard = TRUE)
#'
#' plot(result0, level = 0.90)
#'
#' data(CreditMDR, package = "DiSSMod")
#' n1 <- 600
#' set.n1 <- sample(1:nrow(CreditMDR), n1)
#' reduce_CreditMDR <- CreditMDR[set.n1,]
#' result1 <- DiSSMod(response = MAJORDRG ~ AGE + INCOME + EXP_INC,
#'                    selection = CARDHLDR ~ AGE + INCOME + OWNRENT + ADEPCNT + SELFEMPL,
#'                    data = reduce_CreditMDR, resp.dist="poi", select.dist = "logis",
#'                    alpha = seq(-0.3, 0.3,length.out = 21), standard = FALSE, verbose = 1)
#'
#' plot(result1)
#'
#' @seealso See also \code{\link{DiSSMod}} and \code{\link[graphics]{plot}}.
#' @export
plot.DiSSMod <- function(x, ...)
{
  plot(x$Nalpha, 2*x$st_loglik[!is.na(x$st_loglik)], type="l", xlab = expression(alpha), ylab=expression("Twice relative log-likelihood of" ~~~ alpha))
  points(x$mle_alpha, min(2*x$st_loglik, na.rm = TRUE), col=2, pch=19)
  ci <- confint.DiSSMod(x, ...)
  abline(h=-qchisq(ci$level, 1), lty=2, col=2)
  rug(c(x$mle_alpha, ci$alpha), col=2)
}

#' German doctor first visits data
#'
#' Data is from Riphahn, Wambach and Million (2003), used for studying longitudinal
#' analysis concerning the usage of the German health insurance system. The original
#' data contain a few years data for patients, but we have only for first year.
#'
#' @usage data(DoctorRWM)
#' @format A data frame with 7293 observations of 26 variables as below;
#' \describe{
#'  \item{ID}{identification number (numeric)}
#'  \item{FEMALE}{female or not (categorical)}
#'  \item{YEAR}{year (categorical)}
#'  \item{AGE}{age (numeric)}
#'  \item{HSAT}{health satisfaction coded 0 (low) to 10 (high) (numeric)}
#'  \item{HANDDUM}{person is handicappe or not (categorical)}
#'  \item{HANDPER}{percentage degree of handicap (numeric)}
#'  \item{HHNINC}{monthly household net income (numeric)}
#'  \item{HHKIDS}{child (ren) below age 16 in household (numeric)}
#'  \item{EDUC}{years of schooling (numeric)}
#'  \item{MARRIED}{person is married or not (categorical)}
#'  \item{HAUPTS}{level of schooling (categorical)}
#'  \item{REALS}{level of schooling (categorical)}
#'  \item{FACHHS}{level of schooling (categorical)}
#'  \item{ABITUR}{level of schooling (categorical)}
#'  \item{UNIV}{level of schooling (categorical)}
#'  \item{WORKING}{employed or not (categorical)}
#'  \item{BLUEC}{person is blue collar worker or not (categorical)}
#'  \item{WHITEC}{person is white collar worker or not (categorical)}
#'  \item{SELF}{person is self-employed or not (categorical)}
#'  \item{BEAMT}{civil servant or not (categorical)}
#'  \item{DOCVIS}{number of doctor visits in last 3 months (numeric)}
#'  \item{HOSPVIS}{number of hospital visits last year (numeric)}
#'  \item{PUBLIC}{person is insured in public health insurance or not (categorical)}
#'  \item{ADDON}{person is insured in add-on insurance or not (categorical)}
#'  \item{INCOME_SCALE}{scaled income; original income/1000 (numeric)}
#' }
#' @source Riphahn, R. T., Wambach, A. and Million, A. (2003) Incentive Effects
#' in the Demand for Health Care: A Bivariate Panel Count Data Estimation,
#'  \emph{Journal of Applied Econometrics}, \strong{18}, 4, 387--405.
#'  Published online 8 October 2002.
#' \url{https://doi.org/10.1002/jae.680}
#' @source \url{http://qed.econ.queensu.ca/jae/2003-v18.4/riphahn-wambach-million/}
#' @references Greene, W. H. (2012) \emph{Econometric Analysis, 7th Edition}. Pearson
#' education.
#' @references Azzalini, A., Kim, H.-M. and Kim, H.-J. (2019) Sample selection
#' models for discrete and other non-Gaussian response variables.
#'  \emph{Statistical Methods & Applications}, \strong{28}, 27--56. First online 30 March 2018.
#' \url{https://doi.org/10.1007/s10260-018-0427-1}
"DoctorRWM"

#' Credit cards derogatory reports data
#'
#' Data is originally from Greene (1992), used for studying statistical credit scoring methods.
#'
#' @usage data(CreditMDR)
#' @format A data frame with 13444 observations of 8 variables as below;
#' \describe{
#'  \item{MAJORDRG}{count of major derogatory reports (numeric)}
#'  \item{CARDHLDR}{1 for cardholders, 0 for denied applicants (categorical)}
#'  \item{AGE}{age in years and twelfths of a year (numeric)}
#'  \item{INCOME}{primary income, divided by 10,000 (numeric)}
#'  \item{OWNRENT}{ownRent, individual owns (1) or rents (0) home (categorical)}
#'  \item{ADEPCNT}{not classified yet (numeric)}
#'  \item{SELFEMPL}{self employed; 1=yes, 0=no (categorical)}
#'  \item{EXP_INC}{average expenditure for 12 months/average monthly income (numeric)}
#' }
#' @source Greene, W. H. (1992) \emph{A Statistical Model for Credit
#' Scoring}. Working Paper No. EC-92-29, Department of Economics, Stern School
#' of Business, New York University, 1992.
#' @source \url{http://pages.stern.nyu.edu/~wgreene/Text/Edition7/tablelist8new.htm}
#' @references Greene, W. H. (2012) \emph{Econometric Analysis, 7th Edition}. Pearson
#' education.
#' @references Azzalini, A., Kim, H.-M. and Kim, H.-J. (2019) Sample selection
#' models for discrete and other non-Gaussian response variables.
#'  \emph{Statistical Methods & Applications}, \strong{28}, 27--56. First online 30 March 2018.
#' \url{https://doi.org/10.1007/s10260-018-0427-1}
"CreditMDR"
