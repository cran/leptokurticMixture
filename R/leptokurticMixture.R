



#' Parsimonious model-based clustering with the multivariate elliptical leptokurtic-normal 
#'
#' Performs parsimonious clustering with the multivariate elliptical leptokurtic-normal (MLN). There are 14 possible scale matrix structure and 2 for the kurtosis parameter for a total of 28 models.
#'
#' @param data A n x p matrix of observations.
#' @param G A integer determine the number of components of the mixture model.
#' @param covModels if NULL fit 14 possible scale matrix structures. Otherwise a character vector where each element has length 3. e.g. c("VVV", "EEE") A character of length 4 such as "VVVV", indicating the model; the covariance and beta parameters. The 1st position controls, lambda, the volume; "V" varying across components or "E" equal across components. The 2nd position controls the eigenvalues; V" varying across components, "E" equal across components or "I" the identity matrix. The 3rd  position controls the orientation; "V" varying across components, "E" equal across components or "I" the identity matrix. 
#' @param betaModels set to "V", "E", "B", "F". "V" varying across components, "E" equal across components, "B" consider both "V" & "E", "F" fixed at the maximum value.
#' @param kml a vector of length 3 indicating, the number of k-means starts, number of random starts and the number of EM iterations used for each start 
#' @param label If \code{NULL} then the data has no known groups. If \code{is.integer} then some of the observations have known groups. If \code{label[i]=k} then observation belongs to group  \code{k}. If \code{label[i]=0} then observation has no known group.
#' @param scale.data Should the data be scaled before clustering. The default is TRUE.
#' @param veo "Variables exceed observations". If TRUE, fit the model even though the  number variables in the model exceeds the number of observations. 
#' @param iterMax The maximum number of EM iterations for each model fitted. 
#' @param tol The tol for the stopping rule; lack of progress. The default is 1e-6 but it depends on the data set.
#' @param pprogress If TRUE print the progress of the function.
#' @param method If FP use the fixed point iteration method otherwise if MM use the MM method. 
#' @return A list of
#' \itemize{
#'   \item startobject - A statement on how the models were initialized 
#'   \item gpar - A list of parameter values for the model choosen by the BIC
#'   \item loglik - A vector of the log-likelihoods values 
#'   \item z - A n x G matrix of the posterior probabilities from the model choosen by the BIC
#'   \item map - A vector the maximum a posteriori derived from z
#'   \item BIC - An array with dimensions (G, number of fitted models, 3). The last dimension indices the loglik, number of free parameters and BIC for each fitted model. 
#'   \item bicModel - Information as list on the model choosen by the BIC.
#' }
#' @examples
#' x1 = rmln(n=100, d=4, mu=rep(5,4), diag(4), beta=2)
#' x2 = rmln(n=100, d=4, mu=rep(-5,4), diag(4), beta=2)
#' x = rbind( x1,x2)
#' mlnFit = pmln(data=x, G=2, covModels=c("VVV", "EEE"), betaModels="B")
#' @export
#' @useDynLib leptokurticMixture, .registration = TRUE
#' @importFrom stats cov cov.wt cutree dgamma dist hclust kmeans pgamma qgamma rexp rnorm runif uniroot var
pmln <- function(data=NULL,  G=1:3, covModels=NULL, betaModels="B", kml=c(1,0,1), label=NULL, scale.data=TRUE, veo=FALSE, iterMax=1000, tol=1e-8, pprogress=FALSE, method="FP")  {
  # gpcm
  # pmlmn- parimonous multivariate leperatoic normal 
  
  
  if ( !(method == "FP" | method == "MM") ) stop('method should be either FP or MM')
  if (is.null(data)) stop('Hey, we need some data, please! data is null')
  if (!is.matrix(data)) stop('The data needs to be in matrix form')
  if (!is.numeric(data)) stop('The data is required to be numeric')
  if (nrow(data) == 1) stop('nrow(data) is equal to 1')
  if (ncol(data) == 1) stop('ncol(data) is equal to 1; This function currently only works with multivariate data p > 1')
  if (scale.data) data = scale(data)
  if (any(is.na(data))) stop('No NAs allowed.')
  
  if ( !is.null(iterMax) ) {
      if ( !(is.numeric(iterMax) & (iterMax > 0) ) ) stop("iterMax must be a positive Integer ")
  } else stop("iterMax is null")
  
  if ( !is.null(tol) ) {
    if ( !(is.numeric(tol) & (tol > 0) ) ) stop("tol must be a positive real ")
  } else  stop("tol is null")
  
  if (is.null(G)) stop('G is NULL')
  G = as.integer(ceiling(G))
  if (!is.integer(G)) stop('G is not a integer')
  if ( any(G < 1)) stop('G is not a positive integer')
  
  if (is.null(covModels) )  covModels = c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "EEV", "VEV", "VVV", "EVE", "VVE", "VEE", "EVV")
  if (betaModels == "B")  mnames = c(paste(covModels, "E", sep=""), paste(covModels, "V", sep="") )
  else if (betaModels == "V") mnames = paste(covModels, "V", sep="")
  else if (betaModels == "E") mnames = paste(covModels, "E", sep="")
  else if (betaModels == "F") mnames = paste(covModels, "F", sep="")
  else stop("betaModels needs to be B, V, E or F")  
  
  if (method == "FP") method = 1
  else method = 2
  
  bic = array(0, dim= c(length(G), length(mnames), 3), dimnames=list(G, mnames, c('loglik', "npar", "BIC")) )
  #BIC = matrix(0, nrow=length(G), ncol=length(mnmaes), dimnames=list(G, mnames) )	
  model = NULL; curBIC = Inf; issue=FALSE;
  for (g in 1:length(G)) {
    for (i in 1:length(mnames)) {
      if ( pprogress ) print(c(G[g],mnames[i]))
      if (veo | npar.model(modelname=mnames[i], p=ncol(data), G=G[g]) < nrow(data)) {
        
       # tem = EM(data =temData$x, G = 4, model="VVVV", n=1000, kml=c(10,0,10), epsilon=1e-8)
        #print(mnames[i] )
        #print( tol)
        #print( iterMax )
        
        a =	try( { EM(data=data, G=G[g], n=iterMax, kml=kml, model=mnames[i], label=label, epsilon= tol, estimation = method ) }, silent = TRUE)
        #print(a$loglik)
        #print( !is.nan(a$maxLoglik) )
        #print( !is.na(a$maxLoglik) )
        #print( min( round(diff(a$loglik), 10 )) )
        #print( diff(a$loglik) )
        #print( !is.infinite(a$maxLoglik ) )
        
        if ( length(a) > 1 ) {
          if ( !is.nan(a$maxLoglik) | !is.na(a$maxLoglik) ) { 
            if ( min( round(diff(a$loglik), 10 )) >= 0 & !is.infinite(a$maxLoglik ) ) {
               bic[g,i,1:2] = c(a$maxLoglik, a$numpar  )
               bic[g,i,3] =  -2*bic[g,i,1] + bic[g,i,2]*log(nrow(data)) 
          
              if ( bic[g,i,3] < curBIC) {
                model  = append(a,list(mtype=mnames[i]))
                curBIC = bic[g,i,3] 
              } 
              
            } else issue=TRUE
          } else issue= TRUE
         } else issue =TRUE
        } else issue= TRUE
              
    if (issue) {
      bic[g,i,1:2] = c(NA, npar.model(modelname=mnames[i], p=ncol(data), G=G[g]) )
      bic[g,i,3] =  NA
      issue = FALSE
    }      

      
    }}
  
  startobject= paste( kml[1], " k-means",  kml[2], " random initializations", " with ", kml[3], " emEM iterations", collapse="" )
  
  bicModel = list(G=length(model$gpar$pi), model=model$mtype, bic=curBIC, loglik=model$maxLoglik )
  val = list( startobject= startobject, gpar=model$gpar, loglik=model$loglik, z=model$z, map=model$map, BIC=bic, bicModel= bicModel)
  #	val = list( model=model, BIC=bic, bicModel= bicModel)
  
  class(val) <- "pmln"
  return(val)
}






#' Generate realizations from the multivariate elliptical leptokurtic-normal distribution
#'
#' This function calculates the log cumulative density function for the multivariate-t with scale matrix equal to the identity matrix. It finds the mode and then uses Gaussian quadrature to estimate the integral.
#'
#' @param n number of observations
#' @param d the dimension of the observations
#' @param mu location parameter of length d
#' @param Sigma (d x d) scatter matrix
#' @param beta the concentration parameter
#' @return A (n x d) matrix of realizations
#' @examples
#' x = rmln(n=10, d=4, mu=rep(0,4), diag(4), beta=2)
#' @export
rmln <- function (n = NULL, d = NULL, mu = NULL, Sigma = NULL, beta = NULL) {
  r = sqrt(rlepchisq(n, d, beta))
  z = matrix(rnorm(n * d), n, d)
  z = z/sqrt(rowSums(z^2))
  x = t( (r * z) %*% chol(Sigma) ) + mu
  t(x)
}



#' EM for the finite mixtures of MLN
#'
#' Performs a number of iterations of the EM for the multivariate elliptical leptokurtic-normal (MLN) distribution until the tolerance for the lack progress or the maximum number of iterations is reached. An implementation of parsimonious clustering models via the eigen-decomposition of the scatter matrix and allowing the concentration parameter to be varying, equal or fixed across components.
#'
#' @param data A n x p matrix of observations.
#' @param G A integer determine the number of components of the mixture model.
#' @param model a character of length 4 such as "VVVV", indicating the model; the covariance and beta parameters. The 1st position controls, lambda, the volume; "V" varying across components or "E" equal across components. The 2nd position controls the eigenvalues; V" varying across components, "E" equal across components or "I" the identity matrix. The 3rd  position controls the orientation; "V" varying across components, "E" equal across components or "I" the identity matrix. The 4th position controls the concentration, beta; "V" varying across components, "E" equal across components, "F" fixed at the maximum value.
#' @param kml a vector of length 3 indicating, the number of k-means starts, number of random starts and the number of EM iterations used for each start 
#' @param n The maximum number of EM iterations.
#' @param epsilon The tolerance for the stopping rule; lack of progress. The default is 1e-6 but it depends on the dataset.
#' @param gpar0 A list of model parameters .
#' @param estimation If 1 (default)  use the fixed point iterations and if 2 the MM algorithm.
#' @param label If \code{NULL} then the data has no known groups. If \code{is.integer} then some of the observations have known groups. If \code{label[i]=k} then observation belongs to group  \code{k}. If \code{label[i]=0} then observation has no known group.
#' @return A list with following items
#' \itemize{
#'   \item loglik - A vector of the loglikelihood values
#'   \item gpar - A list containing the parameters values
#'   \item z - A n x G matrix of the posterior probabilities
#'   \item map - A vector the maximum a posteriori derived from z
#'   \item label - The input provided. 
#'   \item numpar - The number of free parameters in the fitted model.
#'   \item maxLoglik - The largest value from loglik.
#' }
#' @examples
#' x1 = rmln(n=100, d=4, mu=rep(5,4), diag(4), beta=2)
#' x2 = rmln(n=100, d=4, mu=rep(-5,4), diag(4), beta=2)
#' x = rbind( x1,x2)
#' mlnFit = EM(data=x, G=2, model="VVVF")
#' @export
#' @useDynLib leptokurticMixture, .registration = TRUE
#' @importFrom stats cov cov.wt cutree dgamma dist hclust kmeans pgamma qgamma rexp rnorm runif uniroot var
EM <- function(data=NULL, G=2, model=NULL, kml=c(1,0,1), n=10, epsilon=1e-2, gpar0=NULL, estimation=1, label=NULL) {
  gpar = gpar0
  if (estimation == 2) {
    if (is.null(gpar)) gpar  = igpar2(data=data, G=G, k=kml[1], m=kml[2], l=kml[3], model=model)
    getlogden   = create.logz2(x=data, G=G, model=model)
    gpar.update = create.update2(x=data, G=G, model=model)
    
  } else {
    if (is.null(gpar)) gpar  = igpar(data=data, G=G, k=kml[1], m=kml[2], l=kml[3], model=model)
    getlogden   = create.logz(x=data, G=G, model=model)
    gpar.update = create.update(x=data, G=G, model=model)
  }
  logden.lab = create.logden.lab(n=nrow(data), G=G, label=label)
  ##

  loglik = numeric(n)
  for (i in 1:3) {
    logden    = getlogden(gpar) + logden.lab
    loglik[i] = logden.to.loglik(logden, gpar$pi)
    wt        = logden.to.weights(logden, gpar$pi)
    gpar      = gpar.update(gpar, wt)
  }

  while ( ( getall(loglik[1:i]) > epsilon) & (i < n ) )  {
    i = i+1
    logden    = getlogden(gpar) + logden.lab
    loglik[i] = logden.to.loglik(logden, gpar$pi)
    wt        = logden.to.weights(logden, gpar$pi)
    gpar      = gpar.update(gpar, wt)
  }
  
  loglik = loglik[2:i]
  logden = getlogden(gpar) + logden.lab
  wt     = logden.to.weights(logden, gpar$pi)
  mapz   = MAP( wt  )
  numpar = npar.model(model, p=ncol(data), G=G)

  val = list(loglik= loglik, gpar=gpar, z=wt, map=mapz, label=label, numpar = numpar, maxLoglik = max(loglik))
  return(val)
}





#' Compare the two methods of estimation
#'
#' Compare the two methods of estimation for fitting a finite mixture of multivariate elliptical leptokurtic-normal distributions; fixed point iterations and MM algorithm.
#'
#' @param mod A character of length 4 such as "VVVV", indicating the model; the covariance and beta parameters.
#' @param data A n x p matrix of observations.
#' @param G The number of components to fit.
#' @param n The maximum number of EM iterations.
#' @param tol The tolerance for the stopping rule; lack of progress. The default is 1e-6 but it depends on the dataset.
#' @param wt a (n x d) matrix of weights for initialization if NULL, then a random weight matrix is generated.
#' @param n0 Given wt, the number of iterations used to obtain the initial parameters
#' @param lab Using given labels (lab) as starting values.
#' @return A vector of times, number of iterations and log-likelihood values.
#'
#' @export
compareEstimation <- function(mod=NULL, data=NULL, G=NULL, n=10^4, tol=1e-6, wt=NULL, n0= 25, lab=NULL) {
  if (is.null(wt)) {
    #wt = matrix(rexp(nrow(data)*G), nrow(data), G)
    #wt  = t(wt/rowSums(wt))
    #wt = lw   = kmeans(x, centers=2, nstart = 20)$cluster
    
    if (is.null(lab)) {
      wt = t(combinewk(weights=matrix(0, nrow=nrow(data), ncol=G), label=kmeans(data, centers=G, nstart = 20, algorithm =  "MacQueen")$cluster))
    } else{
      wt = t(combinewk(weights=matrix(0, nrow=nrow(data), ncol=G), label=lab))
    }
  }

  gpar1 = gpar.wt(data=data, G=G, n=n0, wt=wt, model=mod, estimation = 1 )
  gpar2 = gpar.wt(data=data, G=G, n=n0, wt=wt, model=mod, estimation = 2)

  #gpar1 = gpar.wt(data=data, G=G, n=n0, wt=wt, model=mod, estimation = 1, gpar0 = par2toPar1(gpar2$gpar) )
  #gpar2 = gpar.wt(data=data, G=G, n=n0, wt=wt, model=mod, estimation = 2, gpar0 = par1toPar2(gpar1$gpar) )

  time1 = system.time({ tem1 = EM(data=data, G=G, model=mod, n=n, label =NULL, gpar0=gpar1$gpar, estimation =1, epsilon=tol) }, gcFirst = TRUE)
  time2 = system.time({ tem2 = EM(data=data, G=G, model=mod, n=n, label =NULL, gpar0=gpar2$gpar, estimation =2, epsilon=tol) }, gcFirst = TRUE)

  val = unlist(list(time1=time1[3], time2=time2[3], loglik1= tem1$maxLoglik, loglik2=tem2$maxLoglik,  n1=length(tem1$loglik), n2=length(tem2$loglik) ))
  return(val)
}





###################
# comparison <- function() {}

gpar.wt <- function(data=NULL, G=2, n=10, wt=NULL, model=NULL, estimation=1, gpar0=NULL) {
  if (is.null(wt)) {
    wt = matrix(rexp(nrow(data)*G), nrow(data), G)
    wt  = t(wt/rowSums(wt))
  }
  
  if ( substr(model,3,3) == "I" ) start.mod = "EIIE"
  else  start.mod = "EEEE"
  
  if (estimation == 1) {
    if (is.null(gpar0)) gpar = rgpar(data, G, wt= t(wt), model=model)
    else gpar = gpar0
    mu.update   = create.update.mu(x=data, G=G, model=start.mod)
    gpar.update = create.update(x=data, G=G, model=start.mod, ak=rep(1,G) )
    getlogden   = create.logz(x=data, G=G, model=start.mod )
    
  } else {
    if (is.null(gpar0)) gpar = rgpar2(data, G, wt=t(wt), model=model)
    else gpar = gpar0
    mu.update   = create.update.mu2(x=data, G=G, model=start.mod)
    gpar.update = create.update2(x=data, G=G, model=model)
    getlogden   = create.logz2(x=data, G=G, model=model )
  }
  
  #for (i in 1:n) { 
  #   gpar = mu.update(gpar, wt)
  #   wt     = logden.to.weights(getlogden(gpar), gpar$pi)
  #}  

 for (i in 1:n) {
   gpar = gpar.update(gpar, wt)
   wt     = logden.to.weights(getlogden(gpar), gpar$pi)
 }   
 gpar$model = model  
  
  loglik = logden.to.loglik(getlogden(gpar), gpar$pi)
  return(list(gpar=gpar, loglik=loglik))
}



create.update.mu2 <- function(x=NULL, G=NULL, model=NULL) {
  n = nrow(x)
  d = ncol(x)
  
  if ( substr(model,3,3) == "I" ) {
    malfn <- malIG
    getz  <- getzI
  } else {
    malfn <- malGG2
    getz  <- getzG
  }
  
  function(gpar=NULL, wt=NULL) {
    ng = rowSums(wt)
    gpar$pi = ng/n
    
    M      = abs( sapply(gpar$beta, getM, d=d))
    mal     = t(malfn(x=x, mu=gpar$mu, sig=gpar$sig))
    logfx1 =  (gpar$beta*( g1(mal, d) )/(1 + gpar$beta*g(mal, d)) - 1/2)*wt/M
    
    gpar$mu = t(sapply(1:G, function(k) { colSums( x*logfx1[k,] )/sum( logfx1[k,] ) } ))
    return(gpar)
  }
}


create.update.mu <- function(x=NULL, G=NULL, model=NULL) {
  n = nrow(x)
  d = ncol(x)
  
  if ( substr(model,3,3) == "I" ) malfn   <- malIG
  else malfn   <- malGG
  
  function(gpar=NULL, wt=NULL) {
    ng = rowSums(wt)
    gpar$pi = ng/n
    
    mal     = t(malfn(x=x, mu=gpar$mu, sig=gpar$sig))
    wfx1    = (1 - 2* gpar$beta*g1(mal, d)/( 1 + gpar$beta*g(mal, d) ))*wt
    gpar$mu = t(sapply(1:G, function(k) { colSums( x*wfx1[k,] )/sum( wfx1[k,] ) } ))
    
    return(gpar)
  }
}

gpar.wt.OLD <- function(data=NULL, G=2, n=10, wt=NULL, model=NULL, estimation=1, gpar0=NULL) {
  if (is.null(wt)) {
    wt = matrix(rexp(nrow(data)*G), nrow(data), G)
    wt  = t(wt/rowSums(wt))
  }
  
  if ( substr(model,3,3) == "I" ) start.mod = "EIIE"
  else  start.mod = "EEEE"
  
  if (estimation == 1) {
    if (is.null(gpar0)) gpar = rgpar(data, G, wt= t(wt), model=start.mod)
    else gpar = gpar0
    gpar.update = create.update(x=data, G=G, model=start.mod)
    getlogden   = create.logz(x=data, G=G, model=start.mod )
  } else {
    if (is.null(gpar0)) gpar = rgpar2(data, G, wt=t(wt), model=model)
    else gpar = gpar0

    gpar.update = create.update2(x=data, G=G, model=start.mod)
    getlogden   = create.logz2(x=data, G=G, model=start.mod )
  }
  
  for (i in 1:n) { 
    gpar = gpar.update(gpar, wt)
    wt     = logden.to.weights(getlogden(gpar), gpar$pi)
  }  
    
  gpar$model = model
  
  loglik = logden.to.loglik(getlogden(gpar), gpar$pi)
  return(list(gpar=gpar, loglik=loglik))
}





par1toPar2 <- function(gpar1=NULL) {

  gpar2 = gpar1
  model = gpar1$model
  if (substr(model,3,3) != "I" ) {
    gpar2$sig$eta = gpar1$sig$invS
    gpar2$sig$gam = gpar1$sig$invS
    for (k in 1:length(gpar1$pi)) {
      tem = eigen(gpar1$sig$invS[k,,])
      gpar2$sig$eta[k,,] = tem$vectors %*% diag(sqrt(tem$values) ) %*% t(tem$vectors)
      gpar2$sig$gam[k,,] = tem$vectors
    }
  }
  if (substr(model,3,3) == "E" ) gpar2$sig$gam = gpar1$sig$gam

  gpar2$sig$R = gpar1$sig$S
  gpar2$sig$A = gpar1$sig$S
  gpar2$sig$S = NULL
  gpar2$sig$invS = NULL
  gpar2
}


par2toPar1 <- function(gpar2=NULL) {
  gpar1 = gpar2
  model = gpar2$model
  if (substr(model,3,3) != "I" ) {
    gpar1$sig$invS = gpar2$sig$eta
    gpar1$sig$S = gpar2$sig$eta
    for (k in 1:length(gpar2$pi)) {
      gpar1$sig$invS[k,,] = crossprod( gpar2$sig$eta[k,,] )
      gpar1$sig$S[k,,] = solve(gpar1$sig$invS[k,,])
    }
  } else {
    gpar1$sig$S = gpar2$sig$R
    for (k in 1:length(gpar2$pi))  gpar1$sig$S[k,,] = gpar1$sig$R[k,,]
  }
  gpar1$sig$A = NULL
  gpar1$sig$R = NULL
  return(gpar1)
}


###################
# Generating.Values <- function() {}


dlepchisq <- function(x=NULL,d=NULL,beta=0) {
  coef = c( 1 + beta/8, - beta/4, beta/8)
  den = coef[1]*dgamma(x, d/2, 1/2) + coef[2]*dgamma(x, d/2+1, 1/2) + coef[3]*dgamma(x, d/2+2, 1/2)
  return(den)
}

plepchisq <- function(q=NULL,d=NULL,beta=0) {
  coef = c( 1 + beta/8, - beta/4, beta/8)
  den = coef[1]*pgamma(q, d/2, 1/2) + coef[2]*pgamma(q, d/2+1, 1/2) + coef[3]*pgamma(q, d/2+2, 1/2)
  return(den)
}

qlepchisq <- function(p=NULL,d=NULL,beta=0) {
  if (beta > 0) {
    coef = c( 1 + beta/8, - beta/4, beta/8)
    qstar = qgamma(p, d/2 + 0:2, 1/2)
    rmax  = sum(coef[c(1,3)]*qstar[c(1,3)] )
    rmin = 0
    val = uniroot( function(x) { plepchisq(x, d, beta) - p }, lower=rmin, upper=rmax )$root
  } else val = qgamma(p, d/2, 1/2)
  return(val)
}

rlepchisq <- function(n=NULL, d=NULL, beta=0) {
  u = runif(n)
  sapply(u, function(u0) { qlepchisq(u0, d, beta) } )
}




###################
# WRAPPER <- function() {}


g   <- function(r=NULL, d=NULL) {  ( r^2 - 2*(d+2)*r + d*(d+2) )/(8*d*(d+2)) }
g1  <- function(r=NULL, d0=NULL) {  ( 2*r - 2*(d0+2) + 0 )/(8*d0*(d0+2))  }
g2  <- function(r=NULL, d0=NULL) {  ( 2 )/(8*d0*(d0+2))   }

logf <-  function(r=NULL, beta=NULL, d0=NULL) {   log(1+beta*g(r,d0)) -r/2 }
logf1 <-  function(r=NULL, beta=NULL, d0=NULL) {   beta*( g1(r,d0=d0) )/(1+beta*g(r,d0)) - 1/2 }
logf2 <-  function(r=NULL, beta=NULL, d0=NULL) {   beta* g2(r,d0=d0)/(1+beta*g(r,d0)) - ( beta*( g1(r,d0=d0) )/(1+beta*g(r,d0)) )^2 }



logden.to.loglik <- function(logz=NULL, pi=NULL) {
  val = sum(log(  colSums( exp(logz + log(pi)) ) ))
  return(val)
}

logden.to.weights <- function(logden=NULL, pi=NULL) {
  G = length(pi)
  if (G > 1) {
    logdpi = logden + log(pi)
    maxz   = apply(logdpi, 2, max)
    logdpi = t(logdpi) - maxz
    den    = exp(logdpi)
    w      = den/rowSums(den)
  } else w = matrix(1,nrow=ncol(logden), ncol=G)
  return( t(w) )
}

create.logden.lab <- function(n=NULL, G=NULL, label=NULL)	{
  mat = matrix(0, nrow=G, ncol=n)

  # known is a numeric with
  # 0 if unknown group membership
  # 1,2,3,.. for label of known group
  if (!is.null(label)) {
    neg.Inf0 = c(-Inf,0)
    kw     = label !=0
    for (k in 1:G) {
      labk    = numeric(G)
      labk[k] = 1
      sumk    = sum(label == k)
      mat[ , label == k] = neg.Inf0[labk +1]
    }
  }
  return(mat)
}







igpar <- function(data=NULL, G=NULL, k=1, m=10, l=1, model=NULL) {
  ### m is the number of random starts
  ### n is the number of EM iterations
  if (k < 0 ) stop("k must be non-negative")
  if (m < 0 ) stop("m must be non-negative")
  if (l < 0 ) stop("l must be non-negative")
  
  km = k + m
  if (km == 0) stop("k + m must be > 0")
  
  ## current and max
  llik = numeric(2)
  maxLoglik = numeric(km+1)
  getlogden   = create.logz(x=data, G=G, model=model)
  
  gpar2 = igpar.kmeans(data=data, G=G,  n=l, model=model)
  llik[1:2] = logden.to.loglik(getlogden(gpar2), gpar2$pi)
  if (is.na(llik[2])) llik[2] = -Inf
  
  for (i in 1:k) {
    gpar1 = igpar.kmeans(data=data, G=G,  n=l, model=model)
    llik[1] = logden.to.loglik(getlogden(gpar1), gpar1$pi)
    #print(c("kmeans", i, llik))
    if ( !is.na(llik[1]) ) {
      if (llik[1] > llik[2]) { 
        gpar2 = gpar1
        llik[2] = llik[1]
      }  
    }
  }
  
  if (m > 0) {
    for (i in 1:m) {
      gpar1 = igpar.wt(data=data, G=G,  n=l, model=model)
      llik[1] = logden.to.loglik(getlogden(gpar1), gpar1$pi)
      #print(c("random", i, llik))
      if ( !is.na(llik[1]) ) {
        if (llik[1] > llik[2]) {
          gpar2 = gpar1
          llik[2] = llik[1]
        }  
      }
    }
  }
  return(gpar2)
}

igpar2 <- function(data=NULL, G=NULL, k=1, m=10, l=1, model=NULL) {
  ### m is the number of random starts
  ### n is the number of EM iterations
  if (k < 0 ) stop("k must be non-negative")
  if (m < 0 ) stop("m must be non-negative")
  if (l < 0 ) stop("l must be non-negative")
  
  km = k + m
  if (km == 0) stop("k + m must be > 0")
  
  ## current and max
  llik = numeric(2)
  maxLoglik = numeric(km+1)
  getlogden   = create.logz2(x=data, G=G, model=model)
  
  gpar2 = igpar.kmeans2(data=data, G=G,  n=l, model=model)
  llik[1:2] = logden.to.loglik(getlogden(gpar2), gpar2$pi)
  if (is.na(llik[2])) llik[2] = -Inf
  
  for (i in 1:k) {
    gpar1 = igpar.kmeans2(data=data, G=G,  n=l, model=model)
    llik[1] = logden.to.loglik(getlogden(gpar1), gpar1$pi)
    #print(c("kmeans", i, llik))
    if ( !is.na(llik[1]) ) {
      if (llik[1] > llik[2]) { 
        gpar2 = gpar1
        llik[2] = llik[1]
      }  
    }
  }
  
  if (m > 0) {
    for (i in 1:m) {
      gpar1 = igpar.wt2(data=data, G=G,  n=l, model=model)
      llik[1] = logden.to.loglik(getlogden(gpar1), gpar1$pi)
      #print(c("random", i, llik))
      if ( !is.na(llik[1]) ) {
        if (llik[1] > llik[2]) {
          gpar2 = gpar1
          llik[2] = llik[1]
        }  
      }
    }
  }
  return(gpar2)
}



igpar.wt <- function(data=NULL, G=2, n=10, wt=NULL, model=NULL) {
  if (is.null(wt)) {
    wt = matrix(rexp(nrow(data)*G), nrow(data), G)
    wt  =wt/rowSums(wt)
  }

  if ( substr(model,3,3) == "I" ) start.mod = "EIIE"
  else  start.mod = "EEEE"
  
  gpar = rgpar(data, G, wt=wt, model=model)
  gpar.update = create.update(x=data, G=G, model=start.mod, ak=rep(1,G) )
  getlogden   = create.logz(x=data, G=G, model=start.mod )
  
  try({ 
    for (i in 1:n) {
    gpar = gpar.update(gpar, wt) 
    wt     = logden.to.weights(getlogden(gpar), gpar$pi)
    } 
  }, silent=TRUE)
  gpar$model = model  

  return(gpar)
}



igpar.wt2 <- function(data=NULL, G=2, n=10, wt=NULL, model=NULL) {
  if (is.null(wt)) {
    wt = matrix(rexp(nrow(data)*G), nrow(data), G)
    wt  =wt/rowSums(wt)
  }
  
  if ( substr(model,3,3) == "I" ) start.mod = "EIIE"
  else  start.mod = "EEEE"
  
  gpar = rgpar2(data, G, wt=wt, model=model)
  gpar.update = create.update2(x=data, G=G, model=start.mod )
  getlogden   = create.logz2(x=data, G=G, model=start.mod )
  
  try({ 
    for (i in 1:n) {
      gpar = gpar.update(gpar, wt) 
      wt     = logden.to.weights(getlogden(gpar), gpar$pi)
    } 
  }, silent=TRUE)
  gpar$model = model  
  
  return(gpar)
}


combinewk <- function(weights=NULL, label=NULL)	{
  # known is a numeric with
  # 0 if unknown group membership
  # 1,2,3,.. for label of known group
  if (is.null(label)) stop('label is null')
  kw     = label !=0
  for (j in 1:ncol(weights)) weights[kw,j] = (label == j)[kw]
  return(weights)
}

igpar.lab <- function(data=NULL, G=2, n=10, lab=NULL, model=NULL) {
  z = combinewk(weights=matrix(0, nrow=nrow(data), ncol=G), label=lab)
  gpar = igpar.wt(data=data, G=G, n=n, wt=z, model=model)
  return(gpar)
}

igpar.kmeans <- function(data=NULL, G=NULL, n=1, nstart = 1, model=NULL) {
  lw   = kmeans(data, centers=G, nstart = nstart)$cluster
  gpar = igpar.lab(data=data, G=G, n=n, lab=lw, model=model)
  return(gpar)
}

igpar.hclust <-  function(data=NULL, G=NULL, n=50, model=NULL) {
  lw   = cutree(hclust(dist(data), "average"),k=G)
  gpar = igpar.lab(data=data, G=G, n=n, lab=lw, model=model)
  return(gpar)
}


igpar.lab2 <- function(data=NULL, G=2, n=10, lab=NULL, model=NULL) {
  z = combinewk(weights=matrix(0, nrow=nrow(data), ncol=G), label=lab)
  gpar = igpar.wt2(data=data, G=G, n=n, wt=z, model=model)
  return(gpar)
}

igpar.kmeans2 <- function(data=NULL, G=NULL, n=1, nstart = 1, model=NULL) {
  lw   = kmeans(data, centers=G, nstart = nstart)$cluster
  gpar = igpar.lab2(data=data, G=G, n=n, lab=lw, model=model)
  return(gpar)
}

igpar.hclust2 <-  function(data=NULL, G=NULL, n=50, model=NULL) {
  lw   = cutree(hclust(dist(data), "average"),k=G)
  gpar = igpar.lab2(data=data, G=G, n=n, lab=lw, model=model)
  return(gpar)
}


getall <- function(loglik) {
  if (length(loglik) <3) stop("must have at least 3 likelihood values")
  n = length(loglik)
  lm1 = loglik[n]
  lm  = loglik[(n-1)]
  lm_1  = loglik[(n-2)]
  am = (lm1 - lm)/(lm - lm_1)
  lm1.Inf = lm + (lm1 - lm)/(1-am)
  val = lm1.Inf - lm
  if (is.nan(val)) val=0
  if (val < 0) val= 1
  if ( lm1 < lm ) val=0
  return( val )
}



MAP <- function(w) {
  z = apply(w, 2, function(z) { which(z==max(z)) } )
  z = as.numeric(z)
  return( z)
}




npar.model <- function(modelname=NULL, p=NULL, G=NULL) {
  val = numeric(4)
  val[1] = G-1  ## pi
  val[2] = G*p  ## mu

  ## scatter matrix
  val[3] = ncovpar(modelname= substr(modelname, 1, 3) , p=p, G=G)

  ## beta
  if (substr(modelname, 4, 4) == "V") val[4] = G
  else if (substr(modelname, 4, 4) == "E") val[4] = 1
  else val[4] = 0

  val = sum(val)
  return(val)
}



ncovpar <- function(modelname=NULL, p=NULL, G=NULL) {
  if (is.null(p)) stop("p is null")
  if (is.null(G)) stop("G is null")
  if (is.null(modelname)) stop("modelname is null")

  if (modelname == "EII") npar = 1
  else if (modelname == "VII") npar = G
  else if (modelname == "EEI") npar = p
  else if (modelname == "VEI") npar = p + G -1
  else if (modelname == "EVI") npar = p*G - G +1
  else if (modelname == "VVI") npar = p*G
  else if (modelname == "EEE") npar = p*(p+1)/2
  else if (modelname == "EEV") npar = G*p*(p+1)/2 - (G-1)*p
  else if (modelname == "VEV") npar = G*p*(p+1)/2 - (G-1)*(p-1)
  else if (modelname == "VVV") npar = G*p*(p+1)/2
  else if (modelname == "EVE") npar = p*(p+1)/2 + (G-1)*(p-1)
  else if (modelname == "VVE") npar = p*(p+1)/2 + (G-1)*p
  else if (modelname == "VEE") npar = p*(p+1)/2 + (G-1)
  else if (modelname == "EVV") npar = G*p*(p+1)/2 - (G-1)
  else stop("modelname is not correctly defined")

  return(npar)
}















newBetaF <- function(beta=NULL, d=NULL, wt=NULL, ng=NULL, mal=NULL) {
  ### Fixed, Do not estimate this quantity.
  beta
}

newBetaV <- function(beta=NULL, d=NULL, wt=NULL, ng=NULL, mal=NULL) {
  b1 = g(mal, d)/(1+ beta* g(mal, d) )
  beta.new = beta + rowSums(wt*b1)/rowSums(wt*b1^2)

  beta.new[beta.new < 0]  =0
  beta.new[beta.new > 4*d*(d + 2)/(d + 4)] = 4*d*(d + 2)/(d + 4)

  return(beta.new)
}


newBetaE <- function(beta=NULL, d=NULL, wt=NULL, ng=NULL, mal=NULL) {
  b1 = g(mal, d)/(1+ beta* g(mal, d) )
  beta.new = mean(beta) + sum(wt*b1)/sum(wt*b1^2)

  beta.new[beta.new < 0] = 0
  beta.new[beta.new > 4*d*(d + 2)/(d + 4)] = 4*d*(d + 2)/(d + 4)
  beta.new = rep(beta.new, length(ng))
  return(beta.new)
}








whichBetaUpdate <- function(model=NULL) {
  if (model == "V") fn = newBetaV
  else if (model == "F") fn = newBetaF #fixed
  else if (model == "E") fn = newBetaE
  else stop("The model is ", substr(model,4,4) )
  return(fn)
}






GPAR2 <- function() {}





create.logz2 <- function(x=NULL, G=NULL, model=NULL) {
  n = nrow(x)
  d = ncol(x)

  if (substr(model,3,3) == "I") malfn <- malIG
  else malfn <- malGG2

  function(gpar) {
    ## n x G

    mal = t(malfn(x=x, mu=gpar$mu, sig=gpar$sig))

    tmalb     = log(1 + gpar$beta*g( mal, d) ) - mal/2
    tlogden   = tmalb - 1/2*(d*log(2*pi) + gpar$sig$logdet)

    return(tlogden) #transpose of log density
  }
}



rgpar2 <- function(data=NULL, G=NULL, model=NULL, wt=NULL) {
  d = ncol(data)
  n = nrow(data)

  if (is.null(wt)) {
    wt = matrix( rexp(n*G), nrow=n, ncol=G )
    wt = t(apply(wt, 1, function(z) {z/sum(z)} ))
  }

  mu =matrix(0, G, d)
  sig = list(R=array(0, dim=c(G,d,d) ), A=array(0, dim=c(G,d,d) ))
  for (k in 1:G) {
    tem = cov.wt(data, wt=wt[,k], method="ML")
    mu[k,]     = tem$center
    sig$R[k,,] = diag(diag(cov(data) ))
  }

  eigval  = sapply(1:G, function(k) { apply(data, 2, var) })
  eigvec  = sapply(1:G, function(k) { diag(d) } )

  if (substr(model,1,2) == "EI") sig$lam = matrix( mean(eigval), G, ncol=d)
  else if (substr(model,1,2) == "VI") sig$lam = matrix(apply( eigval,2, mean),G,d)
  else if (substr(model,1,2) == "VV") sig$lam = t(eigval)
  else if (substr(model,1,2) == "EE") sig$lam = matrix(apply(eigval,1,mean), G, ncol=d, byrow=TRUE)
  else if (substr(model,1,2) == "EV") sig$lam = matrix(apply(eigval,1,mean), G, ncol=d, byrow=TRUE)
  else if (substr(model,1,2) == "VE") sig$lam = matrix(apply(eigval,1,mean), G, ncol=d, byrow=TRUE)
  else stop("The given model, ", model, " is not valid at position 2")


  if (substr(model,3,3) == "I") sig$gam = NULL
  else if (substr(model,3,3) == "V") sig$gam =  aperm(array( diag(d), dim=c(d,d,G) ), c(3,2,1))
  else if (substr(model,3,3) == "E") sig$gam = eigen(cov.wt(data, method="ML")$cov)$vectors
  else stop("The given model, ", model, " is not valid at position 3")

  sig$logdet = rowSums( log(sig$lam) )

  if (substr(model,3,3) != "I") {
    sig$eta = array(0, dim=c(G,d,d) )
    #for (k in 1:G)  sig$eta[k,,] = solve(sig$R[k,,])
    for (k in 1:G)  sig$eta[k,,] = diag(d)
  }

  if (substr(model,4,4) == "V" ) beta  = rep( 4*d*(d + 2)/(d + 4) , G)# rep( 1, G)
  else if (substr(model,4,4) == "E" ) beta = rep( 4*d*(d + 2)/(d + 4) , G)
  else if (substr(model,4,4) == "F" ) beta = rep( 4*d*(d + 2)/(d + 4) , G)
  else stop("The given model, ", model, " is not valid at position 4")

  val = list(mu=mu, sig=sig, beta=beta, pi=apply(wt,2,mean), model=model )

  return(val)
}



whichSigUpdate2 <- function(model=NULL) {
  if (model == "VII") fn = getVII2
  else if (model == "EII") fn = getEII2
  else if (model == "EVI") fn = getEVI2
  else if (model == "VVI") fn = getVVI2
  else if (model == "EEI") fn = getEEI2
  else if (model == "VVV") fn = getVVV2
  else if (model == "EEE") fn = getEEE2
  else if (model == "VEI") fn = getVEI2
  else if (model == "EVV") fn = getEVV2
  else if (model == "VEV") fn = getVEV2
  else if (model == "EEV") fn = getEEV2
  else if (model == "VEE") fn = getVEE2
  else if (model == "EVE") fn = getEVE2
  else if (model == "VVE") fn = getVVE2
  else stop("The model is ", substr(model,1,3) )
  fn
}




malIG <- function(x=NULL, mu=NULL, sig=NULL) {
  sapply(1:nrow(mu), function(k) {
    xr2 = ( t(x) - mu[k,] )^2/sig$lam[k,]
    colSums( xr2 )
  })
}

malGG2 <- function(x=NULL, mu=NULL, sig=NULL) {
  sapply(1:nrow(mu), function(k) {
    tx.mu =  t(x) - mu[k,]
    colSums( (sig$eta[k,,] %*% tx.mu)^2 )
  })
}


getRoot <- function(beta=NULL, d=NULL) {
  coefz = numeric(5)
  coefz[1] = -3*d*(d + 2)^2*(beta + 8)*(beta*d + 4*beta - 8*d)
  coefz[2] =  8*beta*(d + 2)^2*(beta*d - beta + 12*d)
  coefz[3] = -6*beta*(d + 2)*(beta*d - 2*beta + 16*d)
  coefz[4] = 0
  coefz[5] = beta^2
  rootz = Re(polyroot(coefz))
  rootz
}

getM <- function(beta=NULL, d=NULL) {
  rootz = getRoot(beta, d)
  min(logf2(rootz, beta, d)*4*rootz + logf1(rootz, beta, d)*2)
}

getzI <- function(tx.mu=NULL, lam=NULL, eta=NULL ) {
  t(tx.mu/sqrt(lam))
}

getzG <- function(tx.mu=NULL, lam=NULL, eta=NULL ) {
  t(eta %*% tx.mu )
}

create.update2  <- function(x=NULL, G=NULL, model=NULL) {
  n = nrow(x)
  d = ncol(x)

  if ( substr(model,3,3) == "I" ) {
    malfn <- malIG
    getz  <- getzI
  } else {
    malfn <- malGG2
    getz  <- getzG
  }

  getBeta = whichBetaUpdate(substr(model,4,4))
  getSig  = whichSigUpdate2(substr(model,1,3))

  function(gpar=NULL, wt=NULL) {
    ng = rowSums(wt)
    gpar$pi = ng/n

    M      = abs( sapply(gpar$beta, getM, d=d))
    mal     = t(malfn(x=x, mu=gpar$mu, sig=gpar$sig))
    logfx1 =  (gpar$beta*( g1(mal, d) )/(1 + gpar$beta*g(mal, d)) - 1/2)*wt/M

    gpar$mu = t(sapply(1:G, function(k) { colSums( x*logfx1[k,] )/sum( logfx1[k,] ) } ))

    mal    = t(malfn(x=x, mu=gpar$mu, sig=gpar$sig))
    fx1    = 1- 2*gpar$beta* g1(mal, d)/( 1 + gpar$beta*g(mal, d) )
    logfx1 =  gpar$beta*( g1(mal, d) )/(1 + gpar$beta*g(mal, d)) - 1/2
    Mfx1   =  (M + 2*logfx1 )*wt

    for (k in 1:G) {
      tx.mu = t(x) -  gpar$mu[k,]
      gpar$sig$R[k,,] = M[k]*( tx.mu %*% ( wt[k,]* t(tx.mu) ) )/sum(wt[k,])
      z0              = Mfx1[k,]* getz( tx.mu, gpar$sig$lam[k,], gpar$sig$eta[k,,] )
      gpar$sig$A[k,,] = ( t(z0) %*% ( t(tx.mu) ))/sum(wt[k,])
    }

    gpar$sig = getSig(gpar$sig, gpar$pi)
    mal = t(malfn(x=x, mu=gpar$mu, sig=gpar$sig))
    gpar$beta = getBeta(beta=gpar$beta, d=d, wt=wt, ng=ng, mal=mal)

    return(gpar)
  }
}



getVVV2 <- function(sig=NULL, pi=NULL) {
  ## VVV
  ## output sig, invS and logdet
  d = ncol(sig$eta)
  for (k in 1:length(pi)) {
    sig$eta[k,,]  = updateEta(A= t(sig$A[k,,])/2, Q= diag(d), G=sig$R[k,,])
    #sig$logdet[k] = -2*log(det(sig$eta[k,,]))
    sig$logdet[k] = -log(det( crossprod(sig$eta[k,,]) ))
  }
  sig
}



getEEE2 <- function(sig=NULL, pi=NULL) {
  # EEE

  d = ncol(sig$eta)
  R = colSums(sig$R*pi, 3)
  A = colSums(sig$A*pi, 3)
  Eta = updateEta(A= t(A)/2, Q= diag(d), G=R)

  for (k in 1:length(pi)) sig$eta[k,,] = Eta

  sig$logdet = rep( -1*log(det( crossprod(Eta) )), length(pi) )
  return(sig)
}






getVVI2 <- function(sig=NULL, pi=NULL) {
  # VVI

  diagR = sapply( 1:length(pi), function(k) { diag(sig$R[k,,]) })
  diagA = sapply( 1:length(pi), function(k) { diag(sig$A[k,,]) })

  eta  = negroot( - diagR, diagA, 1)
  sig$lam = t( 1/eta^2 )
  sig$logdet = rowSums(log(sig$lam))
  sig
}



negroot <- function(a=NULL, b=NULL, c=NULL) {
  (-b - sqrt(b^2 - 4*a*c))/(2*a)
}


getVII2 <- function(sig=NULL, pi=NULL) {
  # VII


  d     = ncol(sig$lam)
  diagR = t(sapply( 1:length(pi), function(k) { diag(sig$R[k,,]) }))
  diagA = t(sapply( 1:length(pi), function(k) { diag(sig$A[k,,]) }))

  eta        = negroot( - rowMeans(diagR),  rowMeans(diagA), 1)
  sig$lam    = matrix( 1/eta^2, nrow=length(pi), ncol=d)
  sig$logdet = rowSums(log(sig$lam))

  return(sig)
}


getEEI2 <- function(sig=NULL, pi=NULL) {
  # EEI

  d     = ncol(sig$lam)

  diagR = t(sapply( 1:length(pi), function(k) { diag(sig$R[k,,]) }))
  diagA = t(sapply( 1:length(pi), function(k) { diag(sig$A[k,,]) }))

  eta        = negroot( - colSums(diagR*pi),  colSums(diagA*pi), 1)

  sig$lam    = matrix( 1/eta^2, nrow=length(pi), ncol=d, byrow=TRUE)
  sig$logdet = rowSums(log(sig$lam))
  return(sig)
}

negroot <- function(a=NULL, b=NULL, c=NULL) {
  (-b - sqrt(b^2 - 4*a*c))/(2*a)
}




getEII2 <- function(sig=NULL, pi=NULL) {
  # EII

  d     = ncol(sig$lam)
  diagR = t(sapply( 1:length(pi), function(k) { diag(sig$R[k,,]) }))
  diagA = t(sapply( 1:length(pi), function(k) { diag(sig$A[k,,]) }))

  eta        = negroot( - sum(diagR*pi)/d,  sum(diagA*pi)/d, 1)

  sig$lam    = matrix( 1/eta^2, nrow=length(pi), ncol=d)
  sig$logdet = rowSums(log(sig$lam))
  return(sig)
}







getVEI2 <- function(sig=NULL, pi=NULL) {
  # getVEI

  d      = ncol(sig$lam)
  loglam = rowMeans(log(sig$lam))
  diagR  = t(sapply( 1:length(pi), function(k) { diag(sig$R[k,,]) }))
  diagA  = t(sapply( 1:length(pi), function(k) { diag(sig$A[k,,]) }))

  psi = exp(colMeans(log(sig$lam) - loglam )  )
  lam = exp(loglam)

  psi = newtRap( psi,  -colSums(diagR*pi/lam), colSums(diagA*pi/sqrt(lam)) )
  eta = negroot( - colMeans( t(diagR)/psi ),  colMeans(t(diagA)/sqrt(psi) ) , 1)

  sig$lam    = matrix( psi, nrow=length(pi), ncol=d, byrow = TRUE)/eta^2


  sig$logdet = rowSums(log(sig$lam))

  return(sig)
}


newtRap <- function(psi=NULL, a=NULL, b=NULL) {
  v0    = sum( b/sqrt(psi) + a/2*1/psi )
  x     = -1/2*log(psi[-1])

  gradA = -1*( exp(x)*b[-1] + a[-1]*exp(2*x) ) + exp( - sum(x) )*b[1] + a[1]*exp(2*(-sum(x)) )

  diagA = exp(x)*b[-1] + 2*a[-1]*exp(2*x)
  c2    = exp(  - sum(x) )*b[1] + 2*a[1]*exp(2*(-sum(x) ))
  x     = x  + 0.5*( gradA/diagA  - (1/diagA) * sum(gradA/diagA)/(1/c2+sum(1/diagA) ) )

  psi.new = exp(-2*c(x, -sum(x)) )
  v1    = sum( b/sqrt(psi.new) + a/2*1/psi.new )

  if (v1 > v0) psi.new = psi
  return( psi.new )
}


getEVI2 <- function(sig=NULL, pi=NULL) {
  # getEVI2

  d      = ncol(sig$lam)
  loglam = mean(log(sig$lam))
  diagR  = sapply( 1:length(pi), function(k) { diag(sig$R[k,,]) })
  diagA  = sapply( 1:length(pi), function(k) { diag(sig$A[k,,]) })

  psi = t(exp(log(sig$lam) - loglam))
  lam = exp(loglam)

  eta = negroot( - sum(colMeans( diagR/psi)*pi ),  sum(colMeans( diagA/sqrt(psi) )*pi), 1)
  sig$lam    = matrix( psi, nrow=length(pi), ncol=d, byrow = TRUE)/eta^2

  psi  = sapply(1:length(pi), function(k) {
    newtRap( psi[,k],  - diagR[,k]/lam, diagA[,k]/sqrt(lam) )
  })

  eta = negroot( - sum(colMeans( diagR/psi)*pi ),  sum(colMeans( diagA/sqrt(psi) )*pi), 1)
  sig$lam    = matrix( psi, nrow=length(pi), ncol=d, byrow = TRUE)/eta^2

  sig$logdet = rowSums(log(sig$lam))

  return(sig)
}




newtRap <- function(psi=NULL, a=NULL, b=NULL) {
  v0    = sum( b/sqrt(psi) + a/2*1/psi )
  x     = -1/2*log(psi[-1])

  gradA = -1*( exp(x)*b[-1] + a[-1]*exp(2*x) ) + exp( - sum(x) )*b[1] + a[1]*exp(2*(-sum(x)) )

  diagA = exp(x)*b[-1] + 2*a[-1]*exp(2*x)
  c2    = exp(  - sum(x) )*b[1] + 2*a[1]*exp(2*(-sum(x) ))
  x     = x  + 1*( gradA/diagA  - (1/diagA) * sum(gradA/diagA)/(1/c2+sum(1/diagA) ) )

  psi.new = exp(-2*c(x, -sum(x)) )
  v1    = sum( b/sqrt(psi.new) + a/2*1/psi.new )

  if (v1 < v0) psi.new = psi
  return( psi.new )
}

getEVV2 <- function(sig=NULL, pi=NULL) {

  d      = ncol(sig$lam)
  loglam = mean(log(sig$lam))
  diagR = sapply( 1:length(pi), function(k) { diag( t(sig$gam[k,,]) %*% sig$R[k,,] %*% (sig$gam[k,,]) ) })
  diagA = sapply( 1:length(pi), function(k) { diag( t(sig$gam[k,,]) %*% sig$A[k,,] %*% (sig$gam[k,,]) ) })

  psi = exp(log(sig$lam) - loglam)
  lam = exp(loglam)

  psi  = sapply(1:length(pi), function(k) {
    newtRap( psi[k,],  - diagR[,k]/lam, diagA[,k]/sqrt(lam) )
  })

  eta     = negroot( - sum(colMeans( diagR/psi)*pi ),  sum(colMeans( diagA/sqrt(psi) )*pi), 1)
  sig$lam = matrix( psi, nrow=length(pi), ncol=d, byrow = TRUE)/eta^2

  for (k in 1:length(pi)) {
    sig$gam[k,,] = gamEVV.MM(R=sig$R[k,,], A=sig$A[k,,], lam=sig$lam[k,], D=sig$gam[k,,])
    sig$eta[k,,] = (sig$gam[k,,]) %*% diag(1/sqrt(sig$lam[k,])) %*% t(sig$gam[k,,])
  }

  sig$logdet = rowSums(log(sig$lam))
  return(sig)
}







getVEV2 <- function(sig=NULL, pi=NULL) {
  # VEV

  d      = ncol(sig$lam)
  diagR = t(sapply( 1:length(pi), function(k) { diag( t(sig$gam[k,,]) %*% sig$R[k,,] %*% (sig$gam[k,,]) ) }))
  diagA = t(sapply( 1:length(pi), function(k) { diag( t(sig$gam[k,,]) %*% sig$A[k,,] %*% (sig$gam[k,,]) ) }))

  loglam = rowMeans(log(sig$lam))
  psi = exp(colMeans(log(sig$lam) - loglam )  )
  lam = exp(loglam)

  psi = newtRap( psi,  -colSums(diagR*pi/lam), colSums(diagA*pi/sqrt(lam)) )
  eta = negroot( - colMeans( t(diagR)/psi ),  colMeans(t(diagA)/sqrt(psi) ) , 1)
  sig$lam    = matrix( psi, nrow=length(pi), ncol=d, byrow = TRUE)/eta^2

  for (k in 1:length(pi)) {
    sig$gam[k,,] = gamEVV.MM(R=sig$R[k,,], A=sig$A[k,,], lam=sig$lam[k,], D=sig$gam[k,,])
    sig$eta[k,,] = (sig$gam[k,,]) %*% diag(1/sqrt(sig$lam[k,])) %*% t(sig$gam[k,,])
  }

  sig$logdet = rowSums(log(sig$lam))
  return(sig)
}





gamEVV.MM <- function(R=NULL, A=NULL, lam=NULL, D=NULL) {
  As = (t(A) + A)

  a1 = max(1/lam)
  z1 = - diag(1/lam) %*% t(D) %*%  t(R)  + a1 * t(D) %*% t(R)

  a2 = max(1/sqrt(lam))
  z2 = diag(1/sqrt(lam)) %*% t(D) %*% t(As) + 1*a2 * t(D) %*% t( As )

  tem1 = svd( (z1 + z2)  )
  D    = (tem1$v) %*% t(tem1$u)
  return(D)
}



getEEV2 <- function(sig=NULL, pi=NULL) {
  # EEV

  d     = ncol(sig$lam)
  diagR = t(sapply( 1:length(pi), function(k) { diag( t(sig$gam[k,,]) %*% sig$R[k,,] %*% (sig$gam[k,,]) ) }))
  diagA = t(sapply( 1:length(pi), function(k) { diag( t(sig$gam[k,,]) %*% sig$A[k,,] %*% (sig$gam[k,,]) ) }))

  eta        = negroot( - colSums(diagR*pi),  colSums(diagA*pi), 1)
  sig$lam    = matrix( 1/eta^2, nrow=length(pi), ncol=d, byrow=TRUE)
  sig$logdet = rowSums(log(sig$lam))

  for (k in 1:length(pi)) {
    sig$gam[k,,] = gamEVV.MM(R=sig$R[k,,], A=sig$A[k,,], lam=sig$lam[k,], D=sig$gam[k,,])
    sig$eta[k,,] = (sig$gam[k,,]) %*% diag(1/sqrt(sig$lam[k,])) %*% t(sig$gam[k,,])
  }

  return(sig)
}








getVEE2 <- function(sig=NULL, pi=NULL) {
  # VEE

  d      = ncol(sig$lam)
  loglam = rowMeans(log(sig$lam))
  diagR = t(sapply( 1:length(pi), function(k) { diag( t(sig$gam) %*% sig$R[k,,] %*% (sig$gam) ) }))
  diagA = t(sapply( 1:length(pi), function(k) { diag( t(sig$gam) %*% sig$A[k,,] %*% (sig$gam) ) }))


  psi = exp(colMeans(log(sig$lam) - loglam )  )
  lam = exp(loglam)

  psi = newtRap( psi,  -colSums(diagR*pi/lam), colSums(diagA*pi/sqrt(lam)) )
  eta = negroot( - colMeans( t(diagR)/psi ),  colMeans(t(diagA)/sqrt(psi) ) , 1)
  sig$lam    = matrix( psi, nrow=length(pi), ncol=d, byrow = TRUE)/eta^2

  sig$gam = diag(d)
  sig$gam  = gamVEE.MM(R=sig$R, A=sig$A, lam=sig$lam, pi=pi, D=sig$gam, d=d)

  for (k in 1:length(pi)) {
    sig$eta[k,,] = (sig$gam) %*% diag(1/sqrt(sig$lam[k,])) %*% t(sig$gam)
  }

  sig$logdet = rowSums(log(sig$lam))
  return(sig)
}




getEVE2 <- function(sig=NULL, pi=NULL) {
  # EVE

  d      = ncol(sig$lam)
  loglam = mean(log(sig$lam))
  diagR = sapply( 1:length(pi), function(k) { diag( t(sig$gam) %*% sig$R[k,,] %*% (sig$gam) ) })
  diagA = sapply( 1:length(pi), function(k) { diag( t(sig$gam) %*% sig$A[k,,] %*% (sig$gam) ) })


  psi = exp(log(sig$lam) - loglam)
  lam = exp(loglam)

  psi  = sapply(1:length(pi), function(k) {
    newtRap( psi[k,],  - diagR[,k]/lam, diagA[,k]/sqrt(lam) )
  })

  eta = negroot( - sum(colMeans( diagR/psi)*pi ),  sum(colMeans( diagA/sqrt(psi) )*pi), 1)

  sig$lam    = matrix( psi, nrow=length(pi), ncol=d, byrow = TRUE)/eta^2
  sig$logdet = rowSums(log(sig$lam))

  sig$gam  = gamVEE.MM(R=sig$R, A=sig$A, lam=sig$lam, pi=pi, D=sig$gam, d=d)

  for (k in 1:length(pi)) {
    sig$eta[k,,] = (sig$gam) %*% diag(1/sqrt(sig$lam[k,])) %*% t(sig$gam)
  }

  return(sig)
}





gamVEE.MM <- function(R=NULL, A=NULL, lam=NULL, pi=NULL, D=NULL, d=NULL) {

  m2 = rowSums(sapply(1:length(pi), function(k) {
    a1 = max(1/lam[k,])
    z1 = - diag(1/lam[k,]) %*% t(D) %*% t(R[k,,]*pi[k])  + a1 * t(D) %*% t(R[k,,]*pi[k])

    a2 = max(1/sqrt(lam[k,]) )
    z2 = diag(1/sqrt(lam[k,])) %*% t(D) %*% ( (A[k,,]+t(A[k,,]) )*pi[k])  + 1*a2 * t(D) %*% t((A[k,,]+t(A[k,,]))*pi[k])
    z1 + z2
  }))
  dim(m2) = c(d,d)
  tem2 = svd( m2 )
  D = (tem2$v) %*% t(tem2$u)

  return(D)
}


getVVE2 <- function(sig=NULL, pi=NULL) {
  # VVE

  d     = ncol(sig$lam)
  diagR = sapply( 1:length(pi), function(k) { diag( t(sig$gam) %*% sig$R[k,,] %*% (sig$gam) ) })
  diagA = sapply( 1:length(pi), function(k) { diag( t(sig$gam) %*% sig$A[k,,] %*% (sig$gam) ) })

  eta  = negroot( - diagR, diagA, 1)
  sig$lam = t( 1/eta^2 )

  sig$gam = gamVEE.MM(R=sig$R, A=sig$A, lam=sig$lam, pi=pi, D=sig$gam, d=d)


  for (k in 1:length(pi)) {
    sig$eta[k,,] = (sig$gam) %*% diag(1/sqrt(sig$lam[k,])) %*% t(sig$gam)
  }

  sig$logdet = rowSums(log(sig$lam))
  sig
}







updateEta = function(A=NULL, Q= NULL, G=NULL ) {
  x     <- rbind(cbind(A, -G), cbind(-Q, t(-A)))
  theta <- .Call("riccatiCareSolution", x)
  return(theta)
}


GPAR1 <- function() {}



rgpar <- function(data=NULL, G=NULL, model=NULL, wt=NULL) {
  d = ncol(data)
  n = nrow(data)

  if (is.null(wt)) {
    wt = matrix( rexp(n*G), nrow=n, ncol=G )
    wt = t(apply(wt, 1, function(z) {z/sum(z)} ))
  }

  mu =matrix(0, G, d)
  sig = list(S=array(0, dim=c(G,d,d) ))
  for (k in 1:G) {
    tem = cov.wt(data, wt=wt[,k], method="ML")
    mu[k,]     = tem$center
    sig$S[k,,] = diag(diag(cov(data) ))
  }

  eigval  = sapply(1:G, function(k) { apply(data, 2, var) })
  eigvec  = sapply(1:G, function(k) { diag(d) } )

  if (substr(model,1,2) == "EI") sig$lam = matrix( mean(eigval), G, ncol=d)
  else if (substr(model,1,2) == "VI") sig$lam = matrix(apply( eigval,2, mean),G,d)
  else if (substr(model,1,2) == "VV") sig$lam = t(eigval)
  else if (substr(model,1,2) == "EE") sig$lam = matrix(apply(eigval,1,mean), G, ncol=d, byrow=TRUE)
  else if (substr(model,1,2) == "EV") sig$lam = matrix(apply(eigval,1,mean), G, ncol=d, byrow=TRUE)
  else if (substr(model,1,2) == "VE") sig$lam = matrix(apply(eigval,1,mean), G, ncol=d, byrow=TRUE)
  else stop("The given model, ", model, " is not valid at position 2")


  if (substr(model,3,3) == "I") sig$gam = NULL
  else if (substr(model,3,3) == "V") sig$gam = NULL
  else if (substr(model,3,3) == "E") sig$gam = eigen(cov.wt(data, method="ML")$cov)$vectors
  else stop("The given model, ", model, " is not valid at position 3")

  sig$logdet = rowSums( log(sig$lam) )

  if (substr(model,3,3) != "I") {
    sig$invS = array(0, dim=c(G,d,d) )
    for (k in 1:G)  sig$invS[k,,] = solve(sig$S[k,,])
  }

  if (substr(model,4,4) == "V" ) beta  = rep( 4*d*(d + 2)/(d + 4) , G)# rep( 1, G)
  else if (substr(model,4,4) == "E" ) beta = rep( 4*d*(d + 2)/(d + 4) , G)
  else if (substr(model,4,4) == "F" ) beta = rep( 4*d*(d + 2)/(d + 4) , G)
  else stop("The given model, ", model, " is not valid at position 4")

  val = list(mu=mu, sig=sig, beta=beta, pi=apply(wt,2,mean), model=model )

  return(val)
}



create.logz <- function(x=NULL, G=NULL, model=NULL) {
  n = nrow(x)
  d = ncol(x)

  if (substr(model,3,3) == "I") malfn <- malIG
  else malfn <- malGG

  function(gpar) {
    ## n x G

    mal = t(malfn(x=x, mu=gpar$mu, sig=gpar$sig))

    tmalb     = log(1 + gpar$beta*g( mal, d) ) - mal/2
    tlogden   = tmalb - 1/2*(d*log(2*pi) + gpar$sig$logdet)

    return(tlogden) #transpose of log density
  }
}




# malGG <- function(x=NULL, mu=NULL, sig=NULL) {
#   sapply(1:nrow(mu), function(k) {
#     malG(x=x, mu=mu[k,], lam=sig$lam[k,], invS=sig$invS[k,,])
#   })
# }


malGG <- function(x=NULL, mu=NULL, sig=NULL) {
  sapply(1:nrow(mu), function(k) {
    tx.mu =  t(x) - mu[k,]
    colSums( tx.mu * sig$invS[k,,] %*% tx.mu)
  })
}


whichSigUpdate <- function(model=NULL) {
  if (model == "VII") fn = getVII
  else if (model == "EII") fn = getEII
  else if (model == "EVI") fn = getEVI
  else if (model == "VVI") fn = getVVI
  else if (model == "EEI") fn = getEEI
  else if (model == "VVV") fn = getVVV
  else if (model == "EEE") fn = getEEE
  else if (model == "VEI") fn = getVEI
  else if (model == "EVV") fn = getEVV
  else if (model == "VEV") fn = getVEV
  else if (model == "EEV") fn = getEEV
  else if (model == "VEE") fn = getVEE
  else if (model == "EVE") fn = getEVE
  else if (model == "VVE") fn = getVVE
  else stop("The model is ", substr(model,1,3) )
  fn
}



create.update  <- function(x=NULL, G=NULL, model=NULL, ak=NULL) {
  n = nrow(x)
  d = ncol(x)

  if ( substr(model,3,3) == "I" ) malfn   <- malIG
  else malfn   <- malGG

  getBeta = whichBetaUpdate(substr(model,4,4))
  getSig  = whichSigUpdate(substr(model,1,3))

  function(gpar=NULL, wt=NULL) {
    ng = rowSums(wt)
    gpar$pi = ng/n

    mal     = t(malfn(x=x, mu=gpar$mu, sig=gpar$sig))
    wfx1    = (1 - 2* gpar$beta*g1(mal, d)/( 1 + gpar$beta*g(mal, d) ))*wt
    gpar$mu = t(sapply(1:G, function(k) { colSums( x*wfx1[k,] )/sum( wfx1[k,] ) } ))

    mal  = t(malfn(x=x, mu=gpar$mu, sig=gpar$sig))
    wfx1 = (1 - 2* gpar$beta*g1(mal, d)/( 1 + gpar$beta*g(mal, d) ))*wt

    if (is.null(ak)) ak   = 1/( 1 +  gpar$beta/(4*d*(d + 2)/(d + 4)) )

    for (k in 1:G) {
      tx.mu = t(x) -  gpar$mu[k,]
      gpar$sig$S[k,,] = ( tx.mu %*% ( wfx1[k,]* t(tx.mu) ) )/sum(wt[k,]) * ak[k] + (1-ak[k])*gpar$sig$S[k,,]
    }

    gpar$sig = getSig(gpar$sig, gpar$pi)

    mal = t(malfn(x=x, mu=gpar$mu, sig=gpar$sig))
    gpar$beta = getBeta(beta=gpar$beta, d=d, wt=wt, ng=ng, mal=mal)

    return(gpar)
  }
}


getVVV <- function(sig=NULL, pi=NULL) {
  ## output sig, invS and logdet
  for (k in 1:length(pi)) {
    sig$invS[k,,] = solve(sig$S[k,,])
    sig$logdet[k] = log(det(sig$S[k,,]))
  }
  sig
}

getVVI <- function(sig=NULL, pi=NULL) {
  # Wk is a list of length G with matrix (p x p)
  # ng is a vector of length G of weights of Wk

  #for (k in 1:length(pi)) sig$lam[k,]   = diag(sig$S[k,,])
  sig$lam    = t(sapply( 1:length(pi), function(k) { diag(sig$S[k,,]) }))
  sig$logdet = rowSums(log(sig$lam))
  sig
}


getEEE <- function(sig=NULL, pi=NULL) {
  # Sk is an array of with dim (p x p x G)
  # ng is a vector of length G of weights of Wk

  S = colSums(sig$S*pi, 3)
  invS = solve(S)
  for (k in 1:length(pi)) sig$invS[k,,] = invS

  sig$logdet = rep( log(det(S)), length(pi) )
  return(sig)
}




getEEI <- function(sig=NULL, pi=NULL) {
  # Wk is a list of length G with matrix (p x p)
  # ng is a vector of length G of weights of Wk

  S = diag( colSums(sig$S*pi, 3) )
  sig$lam = matrix(S, nrow=length(pi), ncol=length(S), byrow=TRUE)
  sig$logdet = rep( sum(log(S)), length(pi) )

  return(sig)
}


getVII <- function(sig=NULL, pi=NULL) {
  # Wk is a list of length G with matrix (p x p)
  # ng is a vector of length G of weights of Wk

  lam     = sapply( 1:length(pi), function(k) { mean(diag(sig$S[k,,])) })
  sig$lam = matrix(lam, nrow=length(pi), ncol=ncol(sig$lam) )
  sig$logdet = rowSums(log(sig$lam))

  return(sig)
}



getEII <- function(sig=NULL, pi=NULL) {
  # Wk is a list of length G with matrix (p x p)
  # ng is a vector of length G of weights of Wk

  S = colSums(sig$S*pi, 3)
  #lam = sum(diag(S))/(sum(pi)*d)
  lam = mean(diag(S))

  sig$lam    = matrix(lam, nrow=length(pi), ncol=ncol(sig$lam) )
  sig$logdet = rowSums(log(sig$lam))

  return(sig)
}


getEVI <- function(sig=NULL, pi=NULL) {
  # Wk is a list of length G with matrix (p x p)
  # ng is a vector of length G of weights of Wk

  S = sapply( 1:length(pi), function(k) { diag(sig$S[k,,]*pi[k]) })

  lam = exp( colMeans(log(S)) )
  lambda = t(S)/lam
  lam = sum(lam)

  sig$lam    = lambda*lam
  sig$logdet = rowSums(log(sig$lam))

  return(sig)
}



getVEI <- function(sig=NULL, pi=NULL) {
  # Wk is a list of length G with matrix (p x p)
  # ng is a vector of length G of weights of Wk

  lam = exp(rowMeans(log(sig$lam)))

  W = diag(colSums(sig$S*pi/lam, 3) )
  B = exp( log(W) - mean(log(W)) )

  lam  =  sapply( 1:length(pi), function(k) { mean( diag(sig$S[k,,])/B) })

  sig$lam    = matrix(B, nrow=length(pi), ncol=ncol(sig$lam), byrow = TRUE )*lam
  sig$logdet = rowSums(log(sig$lam))

  return(sig)
}



getEVV <- function(sig=NULL, pi=NULL ) {
  # eplison=1e-12, max.iter= 100;
  # Wk is a list of length G with matrix (p x p)
  # ng is a vector of length G of weights of Wk
  # lam  =  sapply( 1:length(pi), function(k) { mean( diag(sig$S[k,,])/B) })

  d = dim(sig$S)[2]
  G = length(pi)
  lamG  = sapply( 1:G, function(k) { det(sig$S[k,,]*pi[k])^(1/d) })
  lam   = sum(lamG)

  for (k in 1:G) sig$invS[k,,] = 1/lam * solve(sig$S[k,,]*pi[k]/lamG[k] )
  sig$logdet = rep( d*log(lam), G)

  sig
}






getVEV <- function(sig=NULL, pi=NULL) {
  # eplison=1e-14, max.iter= 100
  # Wk is a list of length G with matrix (p x p)
  # ng is a vector of length G of weights of Wk

  val = sig$lam
  d = ncol(sig$lam)
  G = length(pi)

  for (k in 1:length(pi)) {
    tem = eigen(sig$S[k,,] * pi[k], symmetric=TRUE )
    # eigen vectors place holder
    sig$invS[k,,] = tem$vectors
    val[k,]       = tem$values
  }

  lam = apply(sig$lam, 1, prod)^(1/d)
  #for (i in 1:100 ) {
  A   = colSums(val/lam)/prod(colSums(val/lam))^(1/d)
  lam = colMeans(t(val)/A)/pi
  #}
  sig$lam    = matrix(A, nrow=G, ncol=ncol(sig$lam), byrow=TRUE )*lam

  for (k in 1:length(pi)) sig$invS[k,,] =  crossprod( t(sig$invS[k,,])/sqrt(sig$lam[k,]) )
  sig$logdet = rowSums(log(sig$lam))

  return(sig)
}



getEEV <- function(sig=NULL, pi=NULL) {
  # eplison=1e-14, max.iter= 100
  # Wk is a list of length G with matrix (p x p)
  # ng is a vector of length G of weights of Wk

  val = sig$lam
  d = ncol(sig$lam)
  G = length(pi)

  for (k in 1:length(pi)) {
    tem = eigen(sig$S[k,,] * pi[k], symmetric=TRUE )
    # eigen vectors place holder
    sig$invS[k,,] = tem$vectors
    val[k,]       = tem$values
  }

  lam   = colSums(val)

  sig$lam    = matrix(lam, nrow=G, ncol=ncol(sig$lam), byrow=TRUE )

  for (k in 1:length(pi)) sig$invS[k,,] =  crossprod( t(sig$invS[k,,])/sqrt(sig$lam[k,]) )
  sig$logdet = rowSums(log(sig$lam))

  return(sig)
}


getVEE <- function(sig=NULL, pi=NULL) {
  # Wk is a list of length G with matrix (p x p)
  # ng is a vector of length G of weights of Wk
  d = ncol(sig$lam)
  G = length(pi)

  lam = rowMeans(sig$lam)
  C = colSums(sig$S*pi/lam, 3)
  C = C/det(C)^(1/d)
  invC = solve(C)

  lam = sapply(1:G, function(k) { sum( sig$S[k,,]* invC ) } )/d
  sig$lam    = matrix(lam, nrow=G, ncol=ncol(sig$lam))

  for (k in 1:length(pi)) sig$invS[k,,] = invC/lam[k]
  sig$logdet = rowSums(log(sig$lam))

  return(sig)
}




newD3.MM <- function(S=NULL, d=NULL, pi=NULL, Ak=NULL, D=NULL) {
  z1 = rowSums(sapply(1:length(pi), function(k) {
    alpha = max(eigen(S[k,,]*pi[k])$values)
    (1/Ak[k,]) * ( t(D) %*% (S[k,,]*pi[k]) )  - alpha * ( t(D)/Ak[k,]  )
  }))
  dim(z1) = c(d,d)
  tem1 = svd(z1)
  D = (tem1$v) %*% t(tem1$u)
}

newD4.MM <-  function(S=NULL, d=NULL, pi=NULL, Ak=NULL, D=NULL) {
  z2 = rowSums(sapply(1:length(pi), function(k) {
    alpha = max(1/Ak[k,])
    (S[k,,]*pi[k]) %*% (D %*% diag(1/Ak[k,]) )  - alpha *((S[k,,]*pi[k]) %*% (D) )
  }))
  dim(z2) = c(d,d)
  tem2 = svd(z2)
  D  = t( (tem2$v) %*% t(tem2$u) )
  return(D)
}


getEVE <- function(sig=NULL, pi=NULL) {
  # Wk is a list of length G with matrix (p x p)
  # ng is a vector of length G of weights of Wk
  d = ncol(sig$lam)
  G = length(pi)

  #Ak = apply(Wk,3, function(z,D) { diag( t(D) %*% z %*% (D) ) }, D=D )
  D = sig$gam

  Ak = sapply(1:G, function(k) { diag( t(D) %*% (sig$S[k,,]*pi[k]) %*% (D) ) } )
  lam = exp( colMeans(log(Ak)) )
  Ak = t(Ak)/lam
  lam = sum(lam)

  ### Upate gam
  D = newD3.MM(sig$S, d, pi, Ak, D)
  D = newD4.MM(sig$S, d, pi, Ak, D)

  Ak = sapply(1:G, function(k) { diag( t(D) %*% (sig$S[k,,]*pi[k]) %*% (D) ) } )
  lam = exp( colMeans(log(Ak)) )
  Ak = t(Ak)/lam
  lam = sum(lam)

  sig$lam    = matrix(Ak, nrow=G, ncol=ncol(sig$lam)) *lam
  sig$gam = D

  for (k in 1:length(pi)) sig$invS[k,,] =  crossprod( t(D)/sqrt(sig$lam[k,]) )
  sig$logdet = rowSums(log(sig$lam))

  return(sig)
}





getVVE <- function(sig=NULL, pi=NULL) {
  # Wk is a list of length G with matrix (p x p)
  # ng is a vector of length G of weights of Wk
  d = ncol(sig$lam);  G = length(pi);

  D = sig$gam
  Ak = t(sapply(1:G, function(k) { diag( t(D) %*% (sig$S[k,,]) %*% (D) ) } ))

  ### Upate gam
  D = newD3.MM(sig$S, d, pi, Ak, D)
  D = newD4.MM(sig$S, d, pi, Ak, D)

  sig$lam = t(sapply(1:G, function(k) { diag( t(D) %*% (sig$S[k,,]) %*% (D) ) } ))
  sig$gam = D

  for (k in 1:length(pi)) sig$invS[k,,] =  crossprod( t(D)/sqrt(sig$lam[k,]) )
  sig$logdet = rowSums(log(sig$lam))

  return(sig)
}























