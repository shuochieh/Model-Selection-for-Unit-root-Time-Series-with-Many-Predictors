#### ---- packages ----

library(glmnet)
library(Ohit)
library(doParallel)
library(fGarch)
library(RcppEigen)

#### ---- Data-generating function ----

data.gen = function (n, p, type, lx, ly) {
  
  if (type == 0) {
    T.model0 = c(0, 0, 0, 1, 1, rep(0, ly - 5), rep(1, 5), rep(0, p - 5),
                 rep(0, 5), rep(1, 5), rep(0, p - 10), rep(0, p * (lx - 2)))
    
    X = matrix(0, ncol = p, nrow = n)
    y = rep(0, n)
    eta = rnorm(n)
    
    for (t in 2:n) {
      X[t,] = 0.8 * X[(t - 1),] + 1 * eta[t] + rnorm(p)
    }
    
    b.rlv = c(3, 3.75, 4.5, 5.25, 6, 6.75, 7.5, 8.25, 9, 9.25)
    for (t in 6:n) {
      y[t] = 0.45 * y[t - 4] + 0.45 * y[t - 5] + b.rlv[1:5] %*% X[(t - 1), 1:5] +
        b.rlv[6:10] %*% X[(t - 2), 6:10] + rt(1, df = 6)
    }
  }
  
  if (type == 1) {
    T.model0 = c(1, 1, 1, rep(0, ly - 3), 
                 c(1, 1), rep(0, p - 2),
                 c(1, 1), rep(0, p - 2),
                 c(1, 1), rep(0, p - 2),
                 c(1, 1), rep(0, p - 2),
                 rep(0, p * (lx - 4)))
    
    X = matrix(0, nrow = n, ncol = p)
    x_Spec = garchSpec(model = list("omega" = 1, "alpha" = 0.2, "beta" = 0))
    pi1 = c(garchSim(spec = x_Spec, n = n + 1))
    pi2 = c(garchSim(spec = x_Spec, n = n + 1))
    
    for (j in 1:p) {
      if (j == 1) {
        nu = pi1 + rnorm(n + 1, sd = 1)
        for (t in 1:n) {
          X[t,1] = 0.8 * nu[t + 1] + 0.1 * nu[t]
        }
      } else if (j == 2) {
        nu = pi2 + rnorm(n + 1, sd = 1)
        for (t in 1:n) {
          X[t,2] = 0.2 * nu[t + 1] + 0.6 * nu[t]
        }
      } else {
        if (j %% 2 == 0) {
          nu = pi2 + rnorm(n + 1, sd = 1)
          for (t in 1:n) {
            X[t,j] = 0.2 * nu[t + 1] + 0.6 * nu[t]
          }
        } else {
          nu = pi1 + rnorm(n + 1, sd = 1)
          for (t in 1:n) {
            X[t,j] = 0.8 * nu[t + 1] + 0.1 * nu[t]
          }
        }
      }
    }
    
    #### :: Iteration ::
    #b.rlv = c(-1.85, 1.80, 1.35, -1.10, -0.90, 0.75, 0.65, -0.60)  * 5
    b.rlv = c(-7.62, 6.89, 6.72, -6.18, -5.55, 4.47, 3.77, -3.10)
    y     = rep(0, n)
    
    Spec  = garchSpec(model = list("omega" = 0.05, "alpha" = 0.5, "beta" = 0.1))
    error = c(garchSim(Spec, n = n))
    
    for( i in 5:n ){
      y[i] = 1.6 * y[i - 1] - 0.2 * y[i - 2] - 0.4 * y[i - 3] + 
        b.rlv[1:2] %*% X[(i - 1), 1:2] + 
        b.rlv[3:4] %*% X[(i - 2), 1:2] +
        b.rlv[5:6] %*% X[(i - 3), 1:2] +
        b.rlv[7:8] %*% X[(i - 4), 1:2] + 
        error[i]
    }
  }
  
  if (type == 2) {
    T.model0 = c( 1, 0, 0, 1, 0, 1, rep( 0, ly-6 ), rep( 1, 5 ), rep( 0, p-5 ),
                  rep( 0, 5 ), rep( 1, 5 ), rep( 0, p-10 ), rep( 0, p*(lx-2) ) )
    
    X = matrix(0, ncol = p, nrow = n)
    y = rep(0, n)
    eta = rnorm(n)
    
    for (t in 2:n) {
      X[t,] = 0.8 * X[(t - 1),] + 2 * eta[t] + rnorm(p)
    }
    
    b.rlv = c(3, 3.75, 4.5, 5.25, 6, 6.75, 7.5, 8.25, 9, 9.25)
    for( t in 7:n ){
      y[t] = y[t - 1] + 0.45 * y[t - 4] - 0.45 * y[t - 6] + 
        b.rlv[1:5] %*% X[(t-1),1:5] + b.rlv[6:10] %*% X[(t-2),6:10] + 
        rt(1, df = 6)
    }
  }
  
  if (type == 3) {
    T.model0 = c(1, 1, 1, rep(0, ly - 3), 
                 rep(1, 5), rep(0, p - 5), 
                 rep(0, 5), rep(1, 5), rep(0, p - 10), 
                 rep(0, p * (lx - 2)))
    
    X = matrix(0, nrow = n, ncol = p)
    A = toeplitz(c(0.6^(0:7), rep(0, p - 8)))
    
    eta = c(A %*% rt(n = p, df = 13))
    for (t in 3:n) {
      eta_ = c(A %*% rt(n = p, df = 13))
      X[t,] = 0.1 * X[t - 1,] - 0.7 * X[t - 2,] + eta_ + 0.7 * eta
      eta = eta_
    }
    
    b.rlv = c(0.82, -1.03, 1.92, -2.21, 2.42, -2.57, 3.28, -3.54, 3.72, -3.90)
    y     = rep(0, n)
    
    Spec  = garchSpec(model = list("omega" = 0.05,"alpha" = 0.05,"beta" = 0.9))
    error = c(garchSim(Spec, n = n))
    
    for( i in 4:n ){
      y[i] = (2 * cos(0.1) + 0.3) * y[i - 1] - (1 + 2 * 0.3 * cos(0.1)) * y[i - 2] +
        0.3 * y[i - 3] + b.rlv[1:5] %*% X[(i - 1),1:5] +
        b.rlv[6:10] %*% X[(i - 2),6:10] + error[i]
    }
  }
  
  if (type == 4) {
    T.model0 = c(1, 1, 1, rep(0, ly - 3), c(1, 1, 1, 0, 0), rep(0, p - 5),
                 rep(0, 5), c(1, 0, 0, 1, 1), rep(0, p - 10), rep(0, p * (lx - 2)))
    
    #### X is generated from a VAR (4) model ####
    X = matrix(0, nrow = n, ncol = p)
    
    #### :: Coefficient matrices ::
    A1 = matrix(0.15, 5, 5) 
    A2 = matrix(-0.1, 5, 5)
    Phi1 = as.matrix( bdiag( mget( rep( "A1", (p/5) ) ) ) )
    Phi2 = as.matrix( bdiag( mget( rep( "A2", (p/5) ) ) ) )
    
    #### :: Iteration ::
    for( i in 5:n ){
      X[i,] = Phi1 %*% X[(i-1),] + Phi2 %*% X[(i-2),] + rt( p, df=10 )
    }
    b.rlv = c(0.60, -0.75, 0.90, -0, 0, -1.35, 0, 0, 1.80, -1.85)
    y     = rep(0, n)
    
    Spec  = garchSpec( model=list( "omega"=0.05, "alpha"=0.05, "beta"=0.9 ) )
    error = c( garchSim( Spec, n=n ) )
    
    for( i in 4:n ){
      y[i] = 1.6 * y[i-1] - 0.2 * y[i-2] - 0.4 * y[i-3] + 
        b.rlv[1:5] %*% X[(i-1),1:5] + b.rlv[6:10] %*% X[(i-2),6:10] + error[i]
    }
    
  }
  
  if (type == 5) {
    T.model0 = c( 1, 1, 1, rep( 0, ly-3 ), rep( 1, 5 ), rep( 0, p-5 ), 
                  rep( 0, 5 ), rep( 1, 5 ), rep( 0, p-10 ), rep( 0, p*(lx-2) ) )
    
    #### X is generated from a VAR (4) model ####
    X = matrix( 0, nrow=n, ncol=p )
    
    #### :: Coefficient matrices ::
    A1   = matrix( 0.15, 5, 5 ) ; A2 = matrix( -0.1, 5, 5 )
    Phi1 = as.matrix( bdiag( mget( rep( "A1", (p/5) ) ) ) )
    Phi2 = as.matrix( bdiag( mget( rep( "A2", (p/5) ) ) ) )
    
    #### :: Iteration ::
    for( i in 5:n ){
      X[i,] = Phi1 %*% X[(i-1),] + Phi2 %*% X[(i-2),] + rt( p, df=10 )
    }
    
    b.rlv = c( 0.82, -1.03, 1.92, -2.21, 2.42, -2.57, 3.28, -3.54, 3.72, -3.90 )
    y     = rep( 0, n )
    
    Spec  = garchSpec( model=list( "omega"=0.05, "alpha"=0.05, "beta"=0.9 ) )
    error = c( garchSim( Spec, n=n ) )
    
    for( i in 4:n ){
      y[i] = ( 2*cos(0.1) + 0.3 )*y[i-1] - ( 1 + 2*0.3*cos(0.1) )*y[i-2] + 0.3*y[i-3] + 
        b.rlv[1:5] %*% X[(i-1),1:5] + b.rlv[6:10] %*% X[(i-2),6:10] + error[i]
    }
  }
  
  ylag = embed( y[-n], ly+1 )[,-1]
  yy   = embed( y[-n], ly+1 )[,1]
  C    = embed( X[-n,], ly+1 )[,( ncol(X) + 1 ):( (lx+1)*ncol(X) )]
  G    = cbind( ylag, C )
  
  newdata = c( y[(n-1):(n-ly)], tail( embed( X, ly+1 )[,( ncol(X) + 1 ):( (lx+1)*ncol(X) )], 1 ) )
  
  return( list( "y"=yy, "covariates"=C, "ylag"=ylag, "allpred"=G, "newdata"=newdata, "new_y"=y[n], "j.star"=T.model0 ) )
}

#### ---- Algorithms -----
proj = function (y, X, intercept = F) {
  if( intercept==T ){ X = cbind( 1, X ) }
  model = fastLm( y=y, X=X )
  coef  = model$coefficients
  u     = model$residuals
  
  return( list( "residuals"=u, "coef"=coef ) )
}

arFSR.basic = function (y, X, ylag, Kn) {
  p = ncol( X )
  n = nrow( X )
  
  Q = matrix( 0, nrow=n, ncol=Kn )
  W = cbind( qr.Q( qr(ylag) ), Q )
  
  j.hat = rep( NA, Kn )
  for( k in 1:Kn ){
    r        = X - ( W[,1:( ncol(ylag) + k - 1 )] %*% ( t( W[,1:( ncol(ylag) + k - 1 )] ) %*% X ) )
    r.norms  = sqrt( colSums( r^2 ) )
    score    = t(r) %*% y / r.norms
    if( k>1 ){ score[ j.hat[1:(k-1)] ] = 0 }
    j.hat[k] = which.max( abs( score ) )
    W[,( ncol(ylag) + k )] = r[,j.hat[k]]/r.norms[j.hat[k]]
  }
  
  return( j.hat )
}

HDIC.call = function (y, X, w_n, n_ar, n_x, p = NULL, c = 1) {
  ###
  # Input:
  #    y: an observation vector of the response variable (n by 1)
  #    X: a matrix of all predictor variables (n by (q+p))
  #    w_n: the penalty function of two arguments n and p_n.
  #    n_ar: # of AR variables in the model.
  #    n_x: # of exogneous predictors in the model.
  #    p: total # of exognenous candidate variables.
  #    c: a tuning constant.
  # Output:
  #    HDIC: a scalar.
  ###
  
  n = length( y )
  u = proj( y=y, X=X )$residuals
  
  HDIC = n * log(mean(u^2)) + c * (n_ar + n_x) * w_n(n, p)
  
  return( HDIC )
}

arFSR_threestep = function (y, X, ylag, Kn,  # FSR parameters
                            w_n, c = 1) {     # HDIC parameters
  ###
  # Input: 
  #     y: the observation vector of the response variable
  #     X: (n by p) matrix of exogenous variables
  #     ylag: (n by q) matrix of lagged dependent variables
  #     Kn: number of FSR iterations
  #     w_n: penalty function for HDIC, which takes two arguments (n, p)
  #     c: tuning constant
  # Output:
  #     J_FSR: indices selected by FSR
  #     HDIC: a vector of HDICs of each nested model
  #     J_HDIC: indices survived after HDIC
  #     betahat: coefficiets estimated from the HDIC model
  #     HDIC_BM: HDIC of the benchmark model
  ###
  
  n = length( y )
  p = ncol( X )
  if (is.matrix(ylag)) {
    q = ncol(ylag)
  } else if (is.vector(ylag)) {
    q = 1
  } else {
    stop("arFSR_HDIC: ylag should be either matrix or vector.")
  }
  
  # FSR
  if( missing(Kn) ){ Kn = floor(10 * (n / (p^(max(2 / q, (q + 1) / (2 * q)))))^(0.4)) }
  J_FSR = arFSR.basic(y = y, X = X, ylag = ylag, Kn = Kn)
  
  # Stopping rule by HDIC
  HDIC = rep( NA, Kn )
  for(i in 1:Kn){
    HDIC[i] = HDIC.call(y = y, X = cbind(ylag, X[,J_FSR[1:i]]), w_n = w_n,
                        n_ar = q, n_x = i, p = p, c = c)
  }
  J_HDIC = J_FSR[1:which.min(HDIC)]
  
  # Trim X by HDIC
  HDIC_bm   = min(HDIC)
  HDIC_Trim = rep( NA, length(J_HDIC) )
  
  if (length( J_HDIC ) > 1) {
    for (j in 1:length(J_HDIC)) {
      HDIC_Trim[j] = HDIC.call(y = y, X = cbind(ylag, X[,J_HDIC[-j]]), w_n = w_n,
                               n_ar = q, n_x = length(J_HDIC) - 1, p = p, c = c)
    }
    X_Trim = J_HDIC[which( HDIC_Trim > HDIC_bm )]
  }else{
    X_Trim = J_HDIC
  }
  
  # Output
  betahat = proj(y = y, X = cbind(ylag, X[,sort(X_Trim)]), intercept = FALSE)$coef
  
  return(list("J_FSR" = J_FSR, "HDIC" = HDIC, "J_HDIC" = J_HDIC, "X_Trim" = sort(X_Trim),
              "HDIC_BM" = HDIC_bm, "betahat" = betahat))
}

FHTH = function (y, X, ylag, Kn = NULL, w_n, c = 1, # arFSR_threestep parameters
                 d_n, eta = NULL, exo_ar_index,     # Hard-thresholding parameters
                 newdata = NULL) {                  # if provided, prediction on newdata is returned
  ###
  # Input:
  #     y: the observation vector of the response variable
  #     X: (n by p) matrix of exogenous variables
  #     ylag: (n by q) matrix of lagged dependent variables
  #     Kn: number of FSR iterations
  #     w_n: penalty function for HDIC, which takes two arguments (n, p)
  #     c: tuning constant
  #     d_n: a tuning function that diverges with n; used to determine threshold for AR coefficients.
  #     eta: an estimate on how many moments the error term has. E.g. eta=2 means fourth moment exists.
  #          default = 2.
  #     exo_ar_index: a list of vectors representing the indices of each exogenous predictors.
  #                   E.g. [(1,3),(2,4)] means the first predictor is in the 1st and 3rd columns of X, and the 
  #                        2nd and 4th column of X corresponds to the second predictors.
  #     newdata: new data for which prediction is needed.
  # Output:
  #     J_FSR: indices selected by FSR.
  #     J_Trim: combined indices of the selected predictors.
  #     HDIC: a vector of HDICs of each nested model.
  #     J_HDIC: indices survived after HDIC.
  #     HDIC_BM: HDIC of the benchmark model.
  #     H: threshold for AR coefficients.
  #     betahat: coefficiets estimated from the final model.
  #     yhat: prediction on newdata if newdata is provided.
  ###
  n = length(y)
  if (is.vector(ylag)) {
    q = 1
  } else if (is.matrix(ylag)) {
    q = ncol(ylag)
  } else {
    stop("FHTH: ylag should be vector, matrix, or null")
  }
  if (is.null(eta)) {
    eta = 2
  }
  model = arFSR_threestep(y, X, ylag, Kn, w_n, c)
  X_ind = model$X_Trim
  s0 = length(X_ind)
  s0_bar = 0
  for (lst in exo_ar_index) {
    s0_bar = s0_bar + any(X_ind %in% lst)
  }
  H = max(q^(3/2) * n^(-1/2), min((s0 + q)^(1/2), (s0_bar^(1/2) * q^(1/eta)))) * d_n(n) / sqrt(n)
  beta_ar = model$betahat[1:q]
  AR_T = c(1:q)[which(abs(beta_ar) > H)]
  
  new_ylag = ylag[,AR_T]
  new_X = X[,X_ind]
  betahat = proj(y = y, X = cbind(new_ylag, new_X))$coef
  J_Trim = c(AR_T, X_ind + q)
  
  if (is.null(newdata)) {
    return(list("J_FSR" = model$J_FSR, "X_Trim" = model$X_Trim, "Y_Trim" = AR_T, "J_Trim" = J_Trim,
                "J_HDIC" = model$J_HDIC, "HDIC" = model$HDIC, "HDIC_BM" = model$HDIC_BM, "H" = H,
                "betahat" = betahat))
  } else {
    selected_predictor = newdata[J_Trim]
    yhat = selected_predictor %*% betahat
    
    return(list("J_FSR" = model$J_FSR, "X_Trim" = model$X_Trim, "Y_Trim" = AR_T, "J_Trim" = J_Trim,
                "J_HDIC" = model$J_HDIC, "HDIC" = model$HDIC, "HDIC_BM" = model$HDIC_BM, "H" = H,
                "betahat" = betahat, "yhat" = yhat))
  }
}

OGA.2011 = function (y, X, Kn, HDIC="HDBIC", newdata) { 
  model   = Ohit( y=y, X=X, Kn=Kn, HDIC_Type=HDIC, intercept=FALSE )
  betahat = model$betahat_Trim$coefficients[,1]
  yhat    = c(betahat) %*% c( newdata[ model$J_Trim ] ) 
  
  return( list( "yhat"=yhat, "J_Trim"=model$J_Trim, "betahat"=betahat ) )
}

arOGA = function (y, X, ylag, Kn, w_n, c = 1, newdata) {
  n.ylag = ifelse(is.matrix(ylag), ncol(ylag), 1)
  R      = proj(y = y, X = ylag)$residuals
  
  # Forward selection
  model = OGA(y = R, X = X, Kn = Kn)
  J_OGA = model$J_OGA
  Kn    = model$Kn
  
  OGA_betahat = proj(y = y, X = cbind(ylag, X[,sort(J_OGA)]), intercept = FALSE)$coef
  OGA_yhat = c(newdata[c(1:n.ylag, n.ylag + sort(J_OGA))]) %*% OGA_betahat
  
  # HDIC
  HDIC = rep( 0, Kn )
  for( i in 1:Kn ){
    HDIC[i] = HDIC.call(y = y, X = cbind(ylag, X[,J_OGA[1:i]]), w_n = w_n, 
                        n_ar = n.ylag, n_x = i, p = ncol(X), c = c)
  }
  J_HDIC = J_OGA[ 1:which.min(HDIC) ]
  
  # Trim by HDIC
  HDIC.bm = min( HDIC )
  
  # Trim X
  if( length(J_HDIC) > 1 ){
    HDIC_X = rep( NA, length(J_HDIC) )
    for( j in 1:length(J_HDIC) ){
      HDIC_X[j] = HDIC.call(y = y, X = cbind(ylag, X[,J_HDIC[-j]]), w_n = w_n, 
                            n_ar = n.ylag, n_x = (length(J_HDIC) - 1), 
                            p = ncol(X), c = c)
    }
    X_Trim = J_HDIC[which( HDIC_X > HDIC.bm )]
  }else{
    X_Trim = J_HDIC
  }
  
  # Trim Y
  HDIC.bm = HDIC.call(y = y, X = cbind(ylag, X[,X_Trim]), w_n = w_n, 
                      n_ar = n.ylag, n_x = length(X_Trim), p = ncol(X), c = c)
  if( n.ylag > 1 ){
    HDIC_Y = rep( NA, n.ylag )
    for( i in 1:n.ylag ){
      HDIC_Y[i] = HDIC.call(y = y, X = cbind(ylag[,-i], X[,X_Trim]), w_n = w_n, 
                            n_ar = (n.ylag - 1), n_x = length(X_Trim), 
                            p = ncol(X), c = c)
    }
    Y_Trim = which( HDIC_Y > HDIC.bm )
  }else{
    Y_Trim = 1
  }
  
  
  J_Trim = c( Y_Trim, ( n.ylag + sort(X_Trim) ) )
  
  # Output
  if (length(X_Trim) == 0 && length(Y_Trim) == 0) {
    betahat = rep( 0, n.ylag+ncol(X) )
    yhat = mean(y)
  } else {
    betahat = proj(y = y, X = cbind(ylag[,Y_Trim], X[,sort(X_Trim)]), intercept = FALSE)$coef
    yhat = c(newdata[J_Trim]) %*% betahat
  }
  
  return(list("J_OGA" = J_OGA, "J_HDIC" = J_HDIC, "J_Trim" = J_Trim, 
              "X_Trim" = sort(X_Trim), "Y_Trim" = Y_Trim, "beta_hat" = betahat, 
              "yhat" = yhat, "J_OGA_sorted" = sort(J_OGA), "OGA_yhat" = OGA_yhat))
}

biclasso = function( x, y, penalty.factor, intercept=FALSE, exclude=NULL ){
  if( missing( penalty.factor ) ){
    penalty.factor = rep( 1, ncol(x) )
  }
  
  inf.lasso = glmnet( x=x, y=y, family="gaussian", alpha=1, penalty.factor=penalty.factor, 
                      intercept=intercept, exclude=exclude )
  coef      = as.matrix( coef( inf.lasso ) )
  lambda    = inf.lasso$lambda 
  ncoef     = inf.lasso$df
  
  yhat = cbind( 1, x ) %*% coef
  
  residual = y - yhat 
  mse = colMeans( (residual)^2 ) 
  # sse = colSums( (residual)^2)
  nvar = ncoef + 1
  bic = nrow( x )*log( mse ) + log( nrow(x) )*nvar
  best.model = which( bic == min( bic ) )
  
  return(list("BM" = best.model, "coef" = coef[,best.model], "bic" = min(bic)))
}

lasso.model.selection = function(y, X, type = "lasso", newdata, penalty.factor, intercept=FALSE){
  if( missing( penalty.factor ) ){
    penalty.factor = rep( 1, ncol(X) )
  }
  model    = biclasso( x = X, y = y, penalty.factor=penalty.factor )
  coef     = model$coef
  selected = which( coef[-1] != 0)
  
  if( type == "adalasso" ){
    tau     = seq( from=0.1, to=2, length.out=15 )
    adacoef = NULL
    bics    = NULL
    for( i in 1:15 ){
      penalty           = rep( 0, ncol(X) )
      penalty[selected] = ( abs( (coef[-1])[selected] ) )^(tau[i])
      model    = biclasso( X, y, penalty.factor=(penalty)^(-1), intercept=intercept, exclude=which(coef[-1]==0) )
      adacoef  = cbind( adacoef, model$coef )
      bics     = c( bics, model$bic )
    }
    
    coef = adacoef[,which.min(bics)]
  }
  
  lassobeta = which( coef[-1] != 0 ) 
  yhat = c( 1, newdata ) %*% coef
  
  return( list( "selected" = lassobeta, "yhat" = yhat, "coef" = coef ) )
}

#### Miscellaneous ----

judge.ms = function(selection, true.model, ly){
  if( sum(true.model) == 0 ){
    exact = xe = ye = ss = xss = yss = -1
  }
  else{
    exact = ss = xe = ye = xss = yss = 0
    if( sum( abs( selection - true.model ) ) == 0 ){ exact = 1 }
    if( length( which( selection - true.model < 0) ) == 0 ){ ss = 1 }
    
    if( sum( abs( selection[1:ly] - true.model[1:ly] ) ) == 0 ){ ye = 1 }
    if( length( which( selection[1:ly] - true.model[1:ly] < 0 ) ) == 0  ){ yss = 1 }
    
    if( sum( abs( selection[-c(1:ly)] - true.model[-c(1:ly)] ) ) == 0 ){ xe = 1 }
    if( length( which( selection[-c(1:ly)] - true.model[-c(1:ly)] < 0 ) ) == 0 ){ xss = 1 }
  }
  
  true_pos = sum((selection == true.model)[which(true.model == 1)])
  false_pos = sum((selection != true.model)[which(true.model == 0)])
  
  return(c(exact, xe, ye, ss, xss, yss, true_pos, false_pos))
}

#### Simulation program

nsD.pls = function(N, p, sp, lx){
  ly = floor(2 * (N^(0.25)))
  
  mod_sel  = matrix(0, ncol = 8, nrow = 8)
  
  n = N + ly + 1
  
  ## data generating process
  dta = data.gen(n = n, p = p, type = sp, lx = lx, ly = ly)
  exo_ar_index = vector(mode = "list", length = p)
  for (i in 1:p) {
    exo_ar_index[[i]] = seq(from = i, by = p, length.out = lx)
  }
  
  ## Estimation
  c = 0.5
  
  w_n = function (n, p) {
    return(p^0.5)
  }
  
  d_n = function (n) {
    return(0.5)
  }
  
  model_1 = FHTH(y = as.vector(dta$y), X = dta$covariates, ylag = dta$ylag,
                 Kn = 40, w_n = w_n, c = c, d_n = d_n, 
                 eta = 2, exo_ar_index = exo_ar_index, newdata = dta$newdata)
  
  model_2 = arOGA(y = as.vector(dta$y), X = dta$covariates, ylag = dta$ylag, 
                  Kn = 40, w_n = w_n, c = c, newdata = dta$newdata)
  
  model_3 = OGA.2011( y=as.vector(dta$y), X=dta$allpred, newdata=dta$newdata, Kn=40 )
  
  model_4 = lasso.model.selection( y=as.vector(dta$y), X=dta$allpred, newdata=dta$newdata,
                                   penalty.factor=rep( 1, ly+(p*lx) ) )
  model_5 = lasso.model.selection( y=as.vector(dta$y), X=dta$allpred, newdata=dta$newdata,
                                   penalty.factor=rep( 1, ly+(p*lx) ), type="adalasso" )
  model_6 = lasso.model.selection( y=as.vector(dta$y), X=dta$allpred, newdata=dta$newdata,
                                   penalty.factor=c( rep(0,ly), rep(1,p*lx) ) )
  model_7 = lasso.model.selection( y=as.vector(dta$y), X=dta$allpred, newdata=dta$newdata,
                                   penalty.factor=c( rep(0,ly), rep(1,p*lx) ), type="adalasso" )
  
  ## Model selection result
  
  # FSR+HDIC+Trim
  selected = rep(0, ly + (p * lx))
  selected[model_1$J_Trim] = 1
  mod_sel[1,] = judge.ms( selection=selected, true.model=dta$j.star, ly=ly )
  
  # arOGA+HDIC+Trim
  selected = rep(0, ly + (p * lx))
  selected[model_2$J_Trim] = 1
  mod_sel[2,] = judge.ms( selection=selected, true.model=dta$j.star, ly=ly )
  
  # arOGA
  selected = rep(0, ly + (p * lx))
  selected[c(1:ly, ly + model_2$J_OGA_sorted)] = 1
  mod_sel[8,] = judge.ms(selection = selected, true.model = dta$j.star, ly = ly)
  
  # OGA (2011)
  selected = rep( 0, ly + (p*lx) )
  selected[ model_3$J_Trim ] = 1
  mod_sel[3,] = judge.ms( selection=selected, true.model=dta$j.star, ly=ly )
  
  # plain Lasso
  selected = rep( 0, ly + (p*lx) )
  selected[ model_4$selected ] = 1
  mod_sel[4,] = judge.ms( selection=selected, true.model=dta$j.star, ly=ly )
  
  # plain adaLasso
  selected = rep( 0, ly + (p*lx) )
  selected[ model_5$selected ] = 1
  mod_sel[5,] = judge.ms( selection=selected, true.model=dta$j.star, ly=ly )
  
  # unpenalized Lasso
  selected = rep( 0, ly + (p*lx) )
  selected[ model_6$selected ] = 1
  mod_sel[6,] = judge.ms( selection=selected, true.model=dta$j.star, ly=ly )
  
  # unpenalized adaLasso
  selected = rep( 0, ly + (p*lx) )
  selected[ model_7$selected ] = 1
  mod_sel[7,] = judge.ms( selection=selected, true.model=dta$j.star, ly=ly )
  
  ## Prediction result
  MSPEs = rep( 0, 8 )
  
  MSPEs[1] = ( dta$new_y - model_1$yhat )^2
  MSPEs[2] = ( dta$new_y - model_2$yhat )^2
  MSPEs[3] = ( dta$new_y - model_3$yhat )^2
  MSPEs[4] = ( dta$new_y - model_4$yhat )^2
  MSPEs[5] = ( dta$new_y - model_5$yhat )^2
  MSPEs[6] = ( dta$new_y - model_6$yhat )^2
  MSPEs[7] = ( dta$new_y - model_7$yhat )^2
  MSPEs[8] = ( dta$new_y - model_2$OGA_yhat )^2
  
  ## Output
  output = NULL
  for( i in 1:8 ){
    output = c( output, mod_sel[i,], MSPEs[i] )
  }
  
  return( output )
}

#### Computing ----

Nslots = as.numeric( Sys.getenv("SLURM_CPUS_PER_TASK" ) )
print( sprintf( "%d slots were allocated", Nslots ) )
cl = makeCluster( Nslots )
registerDoParallel( cl )

setwd( "/home/shuang42/High Dim Unit Root" )
array.id = as.numeric( Sys.getenv( "SLURM_ARRAY_TASK_ID" ) )

N   = 200
p   = 100
lx  = 4
sp  = 3

start_time = Sys.time()

output = foreach( ii=1:100, .combine=rbind, .inorder=FALSE,
                  .packages=c( "Ohit", "glmnet", "RcppEigen", "fGarch" ) ) %dopar% {
                    sink()
                    cat( "\n", "Starting iteration ", ii, ". Code has been running for ", difftime( Sys.time(), start_time, units="mins" ),
                         "minutes.", "\n", sep="" )
                    nsD.pls(N = N, p = p, sp = sp, lx = lx)
                  }
stopCluster( cl )

#### Wrap up ----

write.csv( as.matrix( output ),
           paste( "DGP_", sp, "_n_", N, "_p_", p, "_lx_", lx, "_id_", array.id, ".csv", sep="" ) )

quit(save = "no")

