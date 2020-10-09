
install.packages("quantreg")
install.packages("cobs")

library(quantreg)
library(cobs)


###########################################################################################################################

# proposed.nonsmooth:
#       proposed method to estimate conditional quantiles, non-smooth version
#       quantile estimation is at the exact estimated nominal tau's of y|y>0 

# Input: 
#       y: n*1 vector, the observed outcome
#       Xl: n*p matrix, the observed covariates in logistic regression 
#       Xq: n*q matrix, the observed covariates in quantile regression
#       xl: m*p matrix, the new covariates in logistic regression, with which conditional quantile function are estimated
#       xq: m*q matrix, the new covariates in quantile regression, with which conditional quantile function are estimated   
#       taus: k*1 vector, the grid of target tau's of y      
#       delta: constant, better to keep the  default 0.499

# Output:
#       Quantiles: m*k matrix, each row is the estimated quantiles for each new case 

###########################################################################################################################

proposed.nonsmooth <- function(y, Xl, Xq, xl, xq, taus=seq(0, 0.99, by=0.01), delta=0.499){
  
  # logistic regression    
  b = 1*(y > 0)
  logistic = glm(b~Xl,family=binomial)
  gamma = logistic$coef
  
  # estimate the probability of observing a positive y for each new case in xl
  xl.temp = cbind(1, xl)
  p_hat = exp(xl.temp %*% gamma)/( 1+exp(xl.temp %*% gamma) )
  
  # size of data, affect the data-driven interpolation window
  n = length(y)
  
  # estimate conditional quantiles for each new case specifically 
  QF = lapply(1:length(p_hat), function(ii){
    
    # one-to-one mapping of the target tau's of y to the nominal tau's of y|y>0
    p_h = p_hat[ii]
    
    part1 = length( which(taus < (1 - as.vector(p_h))) ) # length of An
    part2 = ( taus[ which(taus >= (1 - as.vector(p_h)) & taus <= (1 - as.vector(p_h) + n^(-delta))) ] - 1 + as.vector(p_h) ) * n^delta # transform in Bn
    part3 = ( taus[ which(taus > (1 - as.vector(p_h) + n^(-delta))) ] - 1 + as.vector(p_h) ) / as.vector(p_h) # transform in Cn
    
    taus.s = c(n^(-delta) / p_h, part3) # esimated nominal tau's of y|y>0
    
    if (any(taus.s > 1)){ # problematic case, the probability of observing a positive y is too small 
      
      quant = c( rep(0, part1) ) # simply 0 for all quantiles
    
    } else { # ordinary case
      
      # quantile regression on the estimated nominal tau's of y|y>0 with positive y only
      qcoefs = rq(y~Xq, tau=taus.s, subset=y>0)$coef
      
      # estimate the conditional quantiles by combining estimates on An, Bn and Cn
      xq.temp = c(1, xq[ii, ])
      fittedy = xq.temp %*% qcoefs
      quant = c( rep(0, part1), fittedy[1]*part2, fittedy[-1] )

    }
  })
  
  return(Quantiles=do.call(rbind, QF))
  
}


###########################################################################################################################

# proposed.smooth:
#       proposed method to estimate conditional quantiles, smooth version
#       quantile estimation is at a fixed grid of nominal tau's of y|y>0
#       approximate the estimated nominal tau's to the prefixed grid
#       match and pick the esimation of quantile

# Input: 
#       y: n*1 vector, the observed outcome
#       Xl: n*p matrix, the observed covariates in logistic regression 
#       Xq: n*q matrix, the observed covariates in quantile regression
#       xl: m*p matrix, the new covariates in logistic regression, with which conditional quantile function are estimated
#       xq: m*q matrix, the new covariates in quantile regression, with which conditional quantile function are estimated   
#       taus: k*1 vector, the grid of target tau's of y  
#       taus.s: a fine grid of vector, the fixed grid of nominal tau's of y|y>0, better to keep the default
#               we need kn (number of taus.s) = o(n^0.5), so up too sample size 10^4, the default is good

# Output:
#       Quantiles: m*k matrix, each row is the estimated quantiles for each new case 

###########################################################################################################################

# COBS constraint, He, X., & Ng, P. (1999)
con = rbind(c(0, 0, 0), # f(0) = 0
            c(0, 0, 0), # f(0) = 0
            c(0, 0, 0)) # f(0) = 0

# approximate the estimated nominal tau's to the prefixed grid of nominal tau's, find the matching location on the prefixed grid
locate.tau = function(tau.target, t=seq(0, 0.99, by=0.01)) {
  loc = 0
  if (tau.target>0)
    loc = which(abs(t-tau.target) == min(abs(t-tau.target)))
  return(loc)
}

proposed.smooth <- function(y, Xl, Xq, xl, xq, taus=seq(0, 0.99, by=0.01), taus.s=seq(0, 0.99, by=0.01)){
  
  # logistic regression    
  b = 1*(y > 0)
  logistic = glm(b~Xl,family=binomial)
  gamma = logistic$coef
  
  # approximate the nominal tau's to the prefixed taus.s, and obtain the locations for all new cases in xl
  xl.temp = cbind(1, xl)
  index_estimated = t(apply(xl.temp, 1, function(z){
    
    # estimate the probability of observing a positive y for each new case in xl
    p_hat = exp(z %*% gamma)/( 1+exp(z %*% gamma) )
    # estimate the exact nominal tau's of y|y>0 
    taustar = pmax( (taus-1+as.vector(p_hat)) / as.vector(p_hat), 0)
    # approximate the nominal tau's to the prefixed taus.s, and obtain the locations
    tau.index = unlist( lapply(taustar,locate.tau) )
  
  }))
  
  # quantile regression on a fixed grid of nominal tau's of y|y>0, taus.s[-1], with positive y only
  qcoef = rq(y~Xq, tau=taus.s[-1], subset=y>0)$coef
  
  # estimate positive quantiles for each new case in xq
  qcoefs = cbind(0, qcoef) # Q_y (0 | x, y>0) = 0
  fittedy = cbind(1, xq) %*% qcoefs

  # obtain a grid of smoothed quantiles given Q_y (0 | x, y>0) = 0
  fittedy_smoothed = t(apply(fittedy, 1, function(z){cobs(x=taus.s, y=z, pointwise=con, nknots = 5)$fitted}))
  
  # estimate the conditional quantiles by combining estimates from the zero and positive part
  len_taus = length(taus)
  len_tauss = length(taus.s)
  
  QF = apply(cbind(fittedy_smoothed, index_estimated), 1, function(z){
    
    # extract the estimate on in positive part
    fit = z[1:len_tauss]
    # extract locations 
    index = z[(len_tauss+1):(len_tauss+len_taus)]
    # combine the zero and positive part 
    c(rep(0, sum(index==0)), fit[index[index>0]])
  
  })
  
  return(Quantiles=t(QF))
  
}


###########################################################################################################################

# AQE
#       proposed inference tool to estimate quantile covariate effect
#       suppose Xj is of interest, wish to estimate the effect of changing it from v to u
#       Delta_tau (Xj; u, v) = E_X(-j) Q_y ( tau | Xj=u, X(-j) ) - Q_y ( tau | Xj=v, X(-j) )

# Input: 
#       y: n*1 vector, the observed outcome
#       Xl: n*p matrix, the observed covariates in logistic regression 
#       Xq: n*q matrix, the observed covariates in quantile regression
#       indexl: location of Xj in Xl
#       indexq: location of Xj in Xq
#       value1: v
#       value2: u
#       taus: k*1 vector, the grid of target tau's of y 
#       method: type of proposed estimation method, non-smooth or smooth

# Output:
#       aqe: k*1 vector, estimated AQE of changing Xj from v to u on the k target tau's of y 

# Note:
#     1. no default for value1 and value2
#     2. no default for indexl and indexq, Xj should be contained in both Xl and Xq, otherwise a logistic regression
#         or a quantile regression is enough
#     3. hypothesis testing can be conducted with bootstrapping on the AQE

###########################################################################################################################

AQE <- function(y, Xl, Xq, indexl, indexq, value1, value2, taus=c(0.1, 0.25, 0.5, 0.75, 0.9), method="non-smooth"){
  
  # size of data, used in the averaging step
  n = length(y)

  # create a matrix of new covariates for logistic regression, 1st half: [Xj=value1, X(-j)], 2nd half: [Xj=value2, X(-j)]
  Xl.temp1 = Xl.temp2 = Xl
  Xl.temp1[, indexl] = value1 
  Xl.temp2[, indexl] = value2
  xl = rbind(Xl.temp1, Xl.temp2)
  
  # create a matrix of new covariates for quantile regression, 1st half: [Xj=value1, X(-j)], 2nd half: [Xj=value2, X(-j)]
  Xq.temp1 = Xq.temp2 = Xq
  Xq.temp1[, indexq] = value1 
  Xq.temp2[, indexq] = value2
  xq = rbind(Xq.temp1, Xq.temp2)

  # estimate quantiles for the new covariates
  if (method=="non-smooth"){
    quant = proposed.nonsmooth(y=y, Xl=Xl, Xq=Xq, xl=xl, xq=xq, taus=taus)
  } else if (method=="smooth"){
    quant = proposed.smooth(y=y, Xl=Xl, Xq=Xq, xl=xl, xq=xq, taus=taus)
  } else{
    return("Please choose non-smooth or smooth")
  }
  
  # AQE estimate = average of Q_y ( tau | Xj=value2, X(-j) ) - Q_y ( tau | Xj=value1, X(-j) )
  return(aqe = apply(quant[(n+1):(2*n), ], 2, mean, na.rm=T) - apply(quant[1:n, ], 2, mean, na.rm=T))  

}



### test example

n = 500     # sample size 

# probability of observing a positive y
p = function(x1,x2,gam0=-1.92,gam1=0.19,gam2=0.02) { # gam0=-7,gam1=0.6,gam2=0.05
  lc = gam0 + gam1*x1 + gam2*x2
  exp(lc)/(1+exp(lc))
}

# quantile coefficient functions
bet0 = function(u){ 100*sqrt(u) }
bet1 = function(u){ u^2 - 4*u }
bet2 = function(u){ sin(u) }

# # quantile coefficient functions estimated from NOMAS
# load("beta_meanpdensity.Rdata")

# new covariates to be evaluated
x1 = c(0, 1)
x2 = qnorm(c(0.1, 0.25, 0.5, 0.75, 0.9), 150, 15)
X = cbind(c(rep(x1[1], 5), rep(x1[2], 5)),  rep(x2, 2))

# simulate the data
x1 = rbinom(n, 1, 0.5)   # e.g. gender
x2 = rnorm(n, 150, 15)   # e.g. blood pressure

b = rbinom(n,1,p(x1,x2))

u = runif(n)
w = bet0(u) + bet1(u)*x1 + bet2(u)*x2

y = b*w

# estimate conditional quantiles
proposed.nonsmooth(y=y, Xl=cbind(x1,x2), Xq=cbind(x1,x2), xl=X, xq=X, taus=c(0.1, 0.25, 0.5, 0.75, 0.9))
proposed.smooth(y=y, Xl=cbind(x1,x2), Xq=cbind(x1,x2), xl=X, xq=X, taus=seq(0.1, 0.9, by=0.1))

# estimate average quantile effect 
AQE(y=y, Xl=cbind(x1,x2), Xq=cbind(x1,x2), indexl=1, indexq=1, value1=0, value2=1)
AQE(y=y, Xl=cbind(x1,x2), Xq=cbind(x1,x2), indexl=1, indexq=1, value1=0, value2=1, method="smooth")


