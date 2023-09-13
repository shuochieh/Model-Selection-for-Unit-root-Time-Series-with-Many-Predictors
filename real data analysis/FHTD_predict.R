setwd("./High Dim Unit Root")
rm(list=ls())

# dta = as.matrix(read.csv("cmort.csv", header = TRUE))
# dta = cbind(dta, dta[,2]^2)
# colnames(dta)[4] = "tempr2"

dta = as.matrix(read.csv("HS_stationary2.csv", header = TRUE))
# dta = dta[1:(nrow(dta) - 32),]

library( glmnet ) 
library( Ohit ) 
library( doParallel )
library( RcppEigen )

################################################################################

proj = function (y, X, intercept = F) {
  if( intercept==T ){ X = cbind( 1, X ) }
  model = fastLm( y=y, X=X )
  coef  = model$coefficients
  u     = model$residuals
  
  return( list( "residuals"=u, "coef"=coef ) )
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
  HDIC = rep(NA, Kn + 1)
  for (i in 0:Kn) {
    if (i == 0) {
      HDIC[i + 1] = HDIC.call(y = y, X = ylag, w_n = w_n, n_ar = q, n_x = 0, p = p, c = c)
    } else {
      HDIC[i + 1] = HDIC.call(y = y, X = cbind(ylag, X[,J_FSR[1:i]]), w_n = w_n,
                              n_ar = q, n_x = i, p = p, c = c)
    }
  }
  if (which.min(HDIC) == 1) {
    J_HDIC = NULL
  } else {
    J_HDIC = J_FSR[1:(which.min(HDIC) - 1)]
  }
  
  # HDIC = rep( NA, Kn )
  # for(i in 1:Kn){
  #   HDIC[i] = HDIC.call(y = y, X = cbind(ylag, X[,J_FSR[1:i]]), w_n = w_n,
  #                       n_ar = q, n_x = i, p = p, c = c)
  # }
  # J_HDIC = J_FSR[1:which.min(HDIC)]
  
  # Trim X by HDIC
  if (!is.null(J_HDIC)) {
    HDIC_bm = min(HDIC)
    HDIC_Trim = rep(NA, length(J_HDIC))
    for (j in 1:length(J_HDIC)) {
      if (length(J_HDIC) == 1) {
        HDIC_Trim = HDIC.call(y = y, X = ylag, w_n = w_n, n_ar = q, n_x = 0,
                              p = p, c = c)
        if (HDIC_Trim > HDIC_bm) {
          X_Trim = J_HDIC
        } else {
          X_Trim = NULL
        }
      } else {
        HDIC_Trim[j] = HDIC.call(y = y, X = cbind(ylag, X[,J_HDIC[-j]]), w_n = w_n,
                                 n_ar = q, n_x = length(J_HDIC) - 1, p = p, c = c)
        if (length(which(HDIC_Trim > HDIC_bm)) == 0) {
          X_Trim = NULL
        } else {
          X_Trim = J_HDIC[which(HDIC_Trim > HDIC_bm)]
        }
      }
    }
  } else {
    HDIC_bm = min(HDIC)
    HDIC_Trim = NULL
    X_Trim = NULL
  }
  # if (length( J_HDIC ) > 1) {
  #   for (j in 1:length(J_HDIC)) {
  #     HDIC_Trim[j] = HDIC.call(y = y, X = cbind(ylag, X[,J_HDIC[-j]]), w_n = w_n,
  #                              n_ar = q, n_x = length(J_HDIC) - 1, p = p, c = c)
  #   }
  #   X_Trim = J_HDIC[which( HDIC_Trim > HDIC_bm )]
  # }else{
  #   X_Trim = J_HDIC
  # }
  
  # Output
  if (is.null(X_Trim)) {
    betahat = proj(y = y, X = ylag, intercept = FALSE)$coef
  } else {
    betahat = proj(y = y, X = cbind(ylag, X[,sort(X_Trim)]), intercept = FALSE)$coef
  }
  
  return(list("J_FSR" = J_FSR, "HDIC" = HDIC, "J_HDIC" = J_HDIC, "X_Trim" = sort(X_Trim),
              "HDIC_BM" = HDIC_bm, "betahat" = betahat))
}

FHTH = function (y, X, ylag, Kn = NULL, w_n, c = 1, # arFSR_threestep parameters
                 d_n, eta = 2, exo_ar_index,        # Hard-thresholding parameters
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
  #     exo_ar_index: a list of vectors representing the start-end index of each exogenous predictors.
  #                   E.g. [(1,3),(4,6)] means the first predictor is in column 1-3 of X, and the 
  #                        4-6th column of X corresponds to the second predictors.
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
  if (length(which(abs(beta_ar) > H)) == 0) {
    # print(paste0("FHTH: hard-thresholding value, ", round(H, 4), 
    #              ", is too large, using only 1 lag instead."))
    AR_T = 1
  } else {
    AR_T = c(1:q)[which(abs(beta_ar) > H)]
  }
  
  new_ylag = ylag[,AR_T]
  if (!is.null(X_ind)) {
    new_X = X[,X_ind]
    betahat = proj(y = y, X = cbind(new_ylag, new_X))$coef
    J_Trim = c(AR_T, X_ind + q)
  } else {
    betahat = proj(y = y, X = new_ylag)$coef
    J_Trim = AR_T
  }
  
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

FHTH.val = function (y, X, ylag, exo_ar_index, 
                     Kn, w_n, d_n, eta = NULL, training_por = 0.75, 
                     grid_size = 10, tune_range = c(0.1, 0.7), 
                     intercept = FALSE, newdata) {
  n = length(y)
  n_train = floor(training_por * n)
  
  y_train = y[1:n_train]
  X_train = X[1:n_train,]
  ylag_train = ylag[1:n_train,]
  
  y_train_orig = y_train
  X_train_orig = X_train
  ylag_train_orig = ylag_train
  
  if (intercept) {
    y_bar = mean(y_train)
    x_bar = colMeans(X_train)
    ylag_bar = colMeans(ylag_train)

    y_train = y_train - y_bar
    X_train = t(t(X_train) - x_bar)
    ylag_train = t(t(ylag_train) - ylag_bar)
  }
  
  y_val = y[(n_train + 1):n]
  X_val = X[(n_train + 1):n,]
  ylag_val = ylag[(n_train + 1):n,]

  # if (intercept) {
  #   y_val = y_val - y_bar
  #   X_val = t(t(X_val) - x_bar)
  #   ylag_val = t(t(ylag_val) - ylag_bar)
  # }
  
  L = tune_range[1]
  U = tune_range[2]
  c = seq(from = L, to = U, length.out = grid_size) 
  d = seq(from = L, to = U, length.out = grid_size) 
  models = vector(mode = "list", length = grid_size * grid_size)
  loss = rep(NA, grid_size * grid_size)
  for (i in 1:grid_size) {
    d_ = function(n, d.prime = d[i]) {
      return(d.prime * d_n(n))
    }
    for (j in 1:grid_size) {
      test_model = FHTH(y = y_train, X = X_train, ylag = ylag_train, Kn = Kn, 
                        w_n = w_n, c = c[j], d_n = d_, eta = eta, 
                        exo_ar_index = exo_ar_index)
      X_Trim = test_model$X_Trim
      Y_Trim = test_model$Y_Trim
      if (intercept) {
        betahat = proj(y = y_train_orig, 
                       X = cbind(ylag_train_orig[,Y_Trim], X_train_orig[,X_Trim]),
                       intercept = TRUE)$coef
        loss[(i - 1) * grid_size + j] = sum((y_val - 
                                             cbind(1, ylag_val[,Y_Trim], X_val[,X_Trim]) %*% betahat )^2)
      } else {
        betahat = test_model$betahat
        loss[(i - 1) * grid_size + j] = sum((y_val - 
                                             cbind(ylag_val[,Y_Trim], X_val[,X_Trim]) %*% betahat )^2)
      }
      models[[(i - 1) * grid_size + j]] = test_model
    }
  }
  
  model = models[[which.min(loss)]]
  X_Trim = model$X_Trim
  Y_Trim = model$Y_Trim
  J_Trim = model$J_Trim
  
  if (!is.null(X_Trim)) {
    if (intercept) {
      betahat = proj(y = y, X = cbind(ylag[,Y_Trim], X[,X_Trim]), intercept = TRUE)$coef
      # betahat = proj(y = y_train_orig, 
      #                X = cbind(ylag_train_orig[,Y_Trim], X_train_orig[,X_Trim]), 
      #                intercept = TRUE)$coef
    } else {
      betahat = proj(y, X = cbind(ylag[,Y_Trim], X[,X_Trim]))$coef
      # betahat = proj(y_train_orig, 
      #                X = cbind(ylag_train_orig[,Y_Trim], X_train_orig[,X_Trim]))$coef
    }
  } else {
    if (intercept) {
      betahat = proj(y, X = cbind(ylag[,Y_Trim]), intercept = TRUE)$coef
      # betahat = proj(y = y_train_orig, 
      #                X = cbind(ylag_train_orig[,Y_Trim]), 
      #                intercept = TRUE)$coef
    } else {
      betahat = proj(y, X = cbind(ylag[,Y_Trim]))$coef
      # betahat = proj(y = y_train_orig, 
      #                X = ylag_train_orig[,Y_Trim])$coef
    }
  }
  if (intercept) {
    yhat = c(1, newdata[J_Trim]) %*% betahat
  } else {
    yhat = c(newdata)[J_Trim] %*% betahat
  }
  i_star = (which.min(loss) - 1) %/% grid_size + 1
  j_star = which.min(loss) %% grid_size
  if (j_star == 0) {
    j_star = grid_size
  }
  c_min = c[j_star]
  d_min = d[i_star]
  
  return( list( "J_FSR" = model$J_FSR, "J_HDIC" = model$J_HDIC, 
                "J_Trim" = model$J_Trim, "X_Trim" = model$X_Trim, 
                "Y_Trim" = model$Y_Trim, "HDIC_BM" = model$HDIC_BM, 
                "H" = model$H, "c_min" = c_min, "d_min" = d_min,
                "betahat" = betahat, "yhat" = yhat ) )
}

OGA.2011 = function (y, X, Kn, HDIC = "HDBIC", intercept = FALSE, newdata) { 
  model   = Ohit(y = y, X = X, Kn = Kn, HDIC_Type = HDIC , intercept = intercept)
  betahat = model$betahat_Trim$coefficients[,1]
  if (intercept) {
    yhat = c(betahat) %*% c(1, newdata[model$J_Trim])
  } else {
    yhat = c(betahat) %*% c(newdata[model$J_Trim]) 
  }

  return(list("yhat" = yhat, "J_Trim" = model$J_Trim, "betahat" = betahat))
}

arOGA = function (y, X, ylag, Kn, w_n, c = 1, newdata) {
  n.ylag = ifelse(is.matrix(ylag), ncol(ylag), 1)
  R = proj(y = y, X = ylag)$residuals
  
  # Forward selection
  model = OGA( y=R, X=X, Kn=Kn )
  J_OGA = model$J_OGA
  Kn    = model$Kn
  
  # HDIC
  HDIC = rep( 0, Kn )
  for( i in 1:Kn ){
    HDIC[i] = HDIC.call(y = y, X = cbind(ylag, X[,J_OGA[1:i]]), w_n = w_n, 
                        n_ar = n.ylag, n_x = i, p = ncol(X), c = c)
  }
  J_HDIC = J_OGA[1:which.min(HDIC)]
  
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
  if( length(X_Trim)==0 && length(Y_Trim)==0 ){
    betahat = rep( 0, n.ylag+ncol(X) )
    yhat    = mean(y)
  }else{
    betahat = proj( y=y, X=cbind( ylag[,Y_Trim], X[,sort(X_Trim)] ), intercept=FALSE )$coef
    yhat    = c( newdata[ J_Trim ] ) %*% betahat
  }
  
  return( list( "J_OGA"=J_OGA, "J_HDIC"=J_HDIC, "J_Trim"=J_Trim, "X_Trim"=sort(X_Trim),
                "Y_Trim"=Y_Trim, "beta_hat"=betahat, "yhat"=yhat ) )
}

arOGA.val = function (y, X, ylag, Kn, eta = NULL, w_n, training_por = 0.75, 
                      grid_size = 10, tune_range = c(0.1, 0.7), 
                      intercept = FALSE, newdata) {
  n = length(y)
  n_train = floor(training_por * n)
  
  y_train    = y[1:n_train]
  X_train    = X[1:n_train,]
  ylag_train = ylag[1:n_train,]
  
  y_train_orig = y_train
  X_train_orig = X_train
  ylag_train_orig = ylag_train
  
  if (intercept) {
    y_bar = mean(y_train)
    x_bar = colMeans(X_train)
    ylag_bar = colMeans(ylag_train)

    y_train = y_train - y_bar
    X_train = t(t(X_train) - x_bar)
    ylag_train = t(t(ylag_train) - ylag_bar)
  }
  
  y_val = y[(n_train+1):n]
  X_val = X[(n_train+1):n,]
  ylag_val = ylag[(n_train+1):n,]
  
  # if (intercept) {
  #   y_val = y_val - y_bar
  #   X_val = t(t(X_val) - x_bar)
  #   ylag_val = t(t(ylag_val) - ylag_bar)
  # }
  
  L = tune_range[1]
  U = tune_range[2]
  c = seq(from = L, to = U, length.out = grid_size)
  loss = rep(NA, grid_size)
  models = vector(mode = "list", length = grid_size)
  for (i in 1:grid_size) {
    test_model = arOGA(y = y_train, X = X_train, ylag = ylag_train, Kn = Kn,
                       w_n = w_n, c = c[i], 
                       newdata = rep(0, ncol(ylag) + ncol(X)))
    X_Trim     = test_model$X_Trim
    Y_Trim     = test_model$Y_Trim
    if (intercept) {
      coef = proj(y = y_train_orig, 
                  X = cbind(ylag_train_orig[,Y_Trim], X_train_orig[,X_Trim]),
                  intercept = TRUE)$coef
      loss[i] = sum((y_val - cbind(1, ylag_val[,Y_Trim], X_val[,X_Trim]) %*% coef)^2)
    } else {
      coef = test_model$beta_hat
      loss[i] = sum((y_val - cbind(ylag_val[,Y_Trim], X_val[,X_Trim]) %*% coef)^2)
    }
    models[[i]]   = test_model
  }
  model = models[[which.min(loss)]]
  X_Trim = model$X_Trim
  Y_Trim = model$Y_Trim
  J_Trim = model$J_Trim
  
  if (!is.null(X_Trim)) {
    if (intercept) {
      betahat = proj(y = y, X = cbind(ylag[,Y_Trim], X[,X_Trim]), intercept = TRUE)$coef
    } else {
      betahat = proj(y, X = cbind(ylag[,Y_Trim], X[,X_Trim]))$coef
    }
  } else {
    if (intercept) {
      betahat = proj(y, X = cbind(ylag[,Y_Trim]), intercept = TRUE)$coef
    } else {
      betahat = proj(y, X = cbind(ylag[,Y_Trim]))$coef
    }
  }
  if (intercept) {
    yhat = c(1, newdata[J_Trim]) %*% betahat
  } else {
    yhat = c(newdata)[J_Trim] %*% betahat
  }
  
  # c_min = c[which.min(loss)]
  # model = arOGA(y = y, X = X, ylag = ylag, Kn = Kn, w_n = w_n, c = c_min, newdata = newdata)
  
  return(list("J_OGA" = model$J_OGA, "J_HDIC" = model$J_HDIC, "J_Trim" = J_Trim, 
              "X_Trim" = X_Trim, "Y_Trim" = Y_Trim, "HDIC_BM" = model$HDIC_BM, 
              "beta_hat" = betahat, "yhat" = yhat))
}

biclasso = function (x, y, penalty.factor, intercept = FALSE, exclude = NULL) {
  if( missing( penalty.factor ) ){
    penalty.factor = rep( 1, ncol(x) )
  }
  
  inf.lasso = glmnet( x=x, y=y, family="gaussian", alpha=1, penalty.factor=penalty.factor, 
                      intercept=intercept, exclude=exclude )
  coef      = as.matrix( coef( inf.lasso ) )
  ncoef     = inf.lasso$df
  
  yhat = cbind( 1, x ) %*% coef
  
  residual   = y - yhat 
  mse        = colMeans( (residual)^2 ) 
  sse        = colSums( (residual)^2 )
  nvar       = ncoef #+ 1
  bic        = nrow( x )*log( mse ) + log( nrow(x) )*nvar
  best.model = which( bic == min( bic ) )
  
  return( list( "BM"=best.model, "coef"=coef[,best.model], "bic"=min(bic) ) )
}

lasso.model.selection = function( y, X, type = "lasso", newdata, penalty.factor, intercept=FALSE ){
  
  if (missing(penalty.factor)) {
    penalty.factor = rep(1, ncol(X))
  }
  
  model    = biclasso( x=X, y=y, penalty.factor=penalty.factor, intercept=intercept )
  coef     = model$coef
  selected = which( coef[-1] != 0 )
  
  if( type == "adalasso" ){
    tau     = seq( from=0.1, to=2, length.out=15 )
    adacoef = NULL
    bics    = NULL
    for( i in 1:15 ){
      penalty           = rep( 0, ncol(X) )
      penalty[selected] = ( abs( (coef[-1])[selected] ) )^(tau[i])
      model    = biclasso( x=X, y=y, penalty.factor=(penalty)^(-1), intercept=intercept, exclude=which(coef[-1]==0) )
      adacoef  = cbind( adacoef, model$coef )
      bics     = c( bics, model$bic )
    }
    coef = adacoef[,which.min(bics)]
  }
  
  yhat     = c( 1, newdata ) %*% coef
  
  return(list("selected" = as.vector(selected), "yhat" = yhat, "coef" = coef))
}

ar.aic = function (y, order.max, intercept = FALSE) {
  ic = rep(NA, order.max)
  temp = embed(y, order.max + 1)
  n = nrow(temp)
  for (i in 1:order.max) {
    model = proj(y = temp[,1], X = temp[,2:(1 + i)], intercept = intercept)
    ic[i] = n * log(mean(model$residuals^2)) + 2 * i
  }
  order_star = which.min(ic)
  model = proj(y = temp[,1], X = temp[,2:(1 + order_star)], intercept = intercept)
  coef = model$coef
  
  if (intercept) {
    yhat = c(1, rev(tail(y, order_star))) %*% coef
    yhat = as.numeric(yhat)
  } else {
    yhat = rev(tail(y, order_star)) %*% coef
    yhat = as.numeric(yhat)
  }
  
  return(list("order" = order_star, "yhat" = yhat, "coef" = coef))
}

rolling.pred = function (macrodata, nprev, i, yar = 6, xar = 6, 
                         eta = 2.25, q0 = 2, delta0 = 0.5, training_por = 0.75,
                         grid_size = 10, tune_range = c(0.1, 0.7),
                         intercept = FALSE) {
  
  data.window = macrodata[(1 + nprev - i):(nrow(macrodata) - i),]
  p = ncol(data.window) - 1 
  
  exo_ar_index = vector(mode = "list", length = p)
  for (j in 1:p) {
    exo_ar_index[[j]] = seq(from = j, by = p, length.out = xar)
  }
  
  # training data
  Y.aux    = embed( data.window[,1], max(yar,xar)+1 )
  X.aux    = embed( data.window[,-1], max(yar,xar)+1 )
  y        = Y.aux[,1]
  ylag     = Y.aux[,2:(yar+1)]
  X        = ( X.aux[,-c(1:p)] )[,1:(p*xar)]
  
  
  # predictors
  ylag.new = tail( Y.aux[,1:yar], 1 )
  XX.new   = tail( X.aux[,1:(xar*p)], 1 )
  newdata  = c( as.vector(ylag.new), as.vector(XX.new) )
  
  w_n = function (n, p, q_ = q0, eta_ = eta, delta_ = delta0) {
    #theta_bar = max(2 / (q_ * eta_), (q_ + 1) / (2 * eta_ * q_))
    #return(p^(2 * theta_bar) * (log(n) / 8)^(1 - delta_))
    #return((p^(1 / eta_)) * log(log(n)))
    return(p^(1 / eta_))
  }
  
  d_n = function (n) {
    return(1)
    # return(log(log(n)))
    # return(1.2 * log(log(n)))
  }
  
  model_FSR     = FHTH.val(y = y, X = X, ylag = ylag, 
                           exo_ar_index = exo_ar_index, Kn = 30, w_n = w_n, 
                           d_n = d_n, eta = eta, training_por = training_por, 
                           grid_size = grid_size, tune_range = tune_range, 
                           intercept = intercept, 
                           newdata = newdata)
  model_OGA     = OGA.2011(y = y, X = cbind(ylag, X), Kn = 30, HDIC = "HDBIC", 
                           intercept = intercept, 
                           newdata = newdata)
  model_arOGA   = arOGA.val(y = y, X = X, ylag = ylag, Kn = 30, eta = eta, 
                            w_n = w_n, training_por = training_por, 
                            grid_size = grid_size, tune_range = tune_range, 
                            intercept = intercept, 
                            newdata = newdata)
  model_lasso   = lasso.model.selection(y = y, X = cbind(ylag, X), type = "lasso", newdata = newdata, 
                                        intercept = intercept, 
                                        penalty.factor = c(rep(1, yar), rep(1, ncol(X))))
  model_alasso  = lasso.model.selection( y=y, X=cbind(ylag,X), type="adalasso", newdata=newdata,
                                         intercept = intercept, 
                                         penalty.factor=c(rep(1, yar), rep(1,ncol(X))))
  model_aralasso = lasso.model.selection(y = y, X = cbind(ylag, X), type = "adalasso", newdata = newdata,
                                         intercept = intercept, 
                                         penalty.factor=c(rep(0, yar), rep(1, ncol(X))))
  
  # Benchmark AR model selected by AIC
  model_ar = ar.aic(as.vector(y), order.max = yar, intercept = intercept)
  yhat = model_ar$yhat
  
  # model_ar     = ar( y, aic=T, order.max=yar, method="ols" )
  # yhat         = predict( model_ar, se=FALSE )
  
  cat("iteration:", nprev - i + 1, "c", model_FSR$c_min, 
      "d", model_FSR$d_min, "FSR", model_FSR$J_Trim,
      "Lasso", model_lasso$selected, "AIC", model_ar$order, "\n")
  
  # Analyze AR variables selection
  fsr_ar_var = rep(0, yar)
  oga_ar_var = rep(0, yar)
  aroga_ar_var = rep(0, yar)
  lasso_ar_var = rep(0, yar)
  alasso_ar_var = rep(0, yar)
  aralasso_ar_var = rep(0, yar)
  araic_ar_var = rep(0, yar)
  
  if (length(model_FSR$Y_Trim) > 0) {
    fsr_ar_var[model_FSR$Y_Trim] = 1
  }
  temp = model_OGA$J_Trim[which(model_OGA$J_Trim <= yar)]
  if (length(temp) > 0) {
    oga_ar_var[temp] = 1
  }
  if (length(model_arOGA$Y_Trim) > 0) {
    aroga_ar_var[model_arOGA$Y_Trim] = 1
  } 
  temp = model_lasso$selected[which(model_lasso$selected <= yar)]
  if (length(temp) > 0) {
    lasso_ar_var[temp] = 1
  }
  temp = model_alasso$selected[which(model_alasso$selected <= yar)]
  if (length(temp) > 0) {
    alasso_ar_var[temp] = 1
  }
  temp = model_aralasso$selected[which(model_aralasso$selected <= yar)]
  if (length(temp) > 0) {
    aralasso_ar_var[temp] = 1
  }
  araic_ar_var[1:model_ar$order] = 1
  

  return(c(model_FSR$yhat, 
           model_OGA$yhat, 
           model_arOGA$yhat,
           model_lasso$yhat, 
           model_alasso$yhat, 
           model_aralasso$yhat, 
           yhat,
           fsr_ar_var,
           oga_ar_var,
           aroga_ar_var,
           lasso_ar_var,
           alasso_ar_var,
           aralasso_ar_var,
           araic_ar_var
           ))
}

# rolling.pc.pred = function (macrodata, nprev, i, yar = 6, xar = 6, 
#                             eta = 2.25, q0 = 2, delta0 = 0.5, 
#                             training_por = 0.75, grid_size = 10, 
#                             tune_range = c(0.1, 0.7), intercept = FALSE) {
#   
#   data.window = macrodata[(1 + nprev - i):(nrow(macrodata) - i),]
#   pca = prcomp(data.window[,-1], retx = TRUE)
#   data.window = as.matrix(cbind(data.window, pca$retx[,1:5]))
#   
#   p = ncol(data.window) - 1 
#   exo_ar_index = vector(mode = "list", length = p)
#   for (j in 1:p) {
#     exo_ar_index[[j]] = seq(from = j, by = p, length.out = xar)
#   }
#   
#   # training data
#   Y.aux    = embed( data.window[,1], max(yar,xar)+1 )
#   X.aux    = embed( data.window[,-1], max(yar,xar)+1 )
#   y        = Y.aux[,1]
#   ylag     = Y.aux[,2:(yar+1)]
#   X        = ( X.aux[,-c(1:p)] )[,1:(p*xar)]
#   
#   
#   # predictors
#   ylag.new = tail( Y.aux[,1:yar], 1 )
#   XX.new   = tail( X.aux[,1:(xar*p)], 1 )
#   newdata  = c( as.vector(ylag.new), as.vector(XX.new) )
#   
#   w_n = function (n, p, q_ = q0, eta_ = eta, delta_ = delta0) {
#     #theta_bar = max(2 / (q_ * eta_), (q_ + 1) / (2 * eta_ * q_))
#     #return(p^(2 * theta_bar) * (log(n) / 8)^(1 - delta_))
#     #return((p^(1 / eta_)) * log(log(n)))
#     return(p^(1 / eta_))
#   }
#   
#   d_n = function (n) {
#     return(1)
#     # return(log(log(n)))
#     # return(1.2 * log(log(n)))
#   }
#   
#   model_FSR     = FHTH.val(y = y, X = X, ylag = ylag, 
#                            exo_ar_index = exo_ar_index, Kn = 40, w_n = w_n, 
#                            d_n = d_n, eta = eta, training_por = training_por, 
#                            grid_size = grid_size, tune_range = tune_range, 
#                            intercept = intercept, newdata = newdata)
#   model_OGA     = OGA.2011(y = y, X = cbind(ylag, X), Kn = 40, HDIC = "HDBIC", 
#                            intercept = intercept, newdata = newdata)
#   model_arOGA   = arOGA.val(y = y, X = X, ylag = ylag, Kn = 40, eta = eta, 
#                             w_n = w_n, training_por = training_por, 
#                             grid_size = grid_size, tune_range = tune_range, 
#                             intercept = intercept, newdata = newdata)
#   model_lasso   = lasso.model.selection(y = y, X = cbind(ylag, X), type = "lasso", newdata = newdata, 
#                                         intercept = TRUE, 
#                                         penalty.factor = c(rep(1, yar), rep(1, ncol(X))))
#   model_alasso  = lasso.model.selection( y=y, X=cbind(ylag,X), type="adalasso", newdata=newdata,
#                                          intercept = TRUE, 
#                                          penalty.factor=c(rep(1, yar), rep(1,ncol(X))))
#   model_aralasso = lasso.model.selection(y = y, X = cbind(ylag, X), type = "adalasso", newdata = newdata,
#                                          intercept = TRUE, 
#                                          penalty.factor=c(rep(0, yar), rep(1, ncol(X))))
#   
#   # Benchmark AR model selected by AIC
#   model_ar     = ar( y, aic=T, order.max=yar, method="ols" )
#   yhat         = predict( model_ar, se=FALSE )
#   
#   cat("iteration:", nprev - i + 1, "c", model_FSR$c_min, 
#       "d", model_FSR$d_min, "FSR", model_FSR$J_Trim,
#       "AR-ALasso", model_aralasso$selected, "AIC", model_ar$order, "\n")
#   
#   return( c( model_FSR$yhat, model_OGA$yhat, 
#              model_arOGA$yhat,
#              model_lasso$yhat, model_alasso$yhat, 
#              model_aralasso$yhat, yhat ) )
# }

################################################################################
# experiment 1

training_por = 0.8
grid_size = 25
tune_range = c(0.1, 0.7)
yar = 18
xar = 18

nprev = 216
intercept = TRUE

################################################################################

no_cluster = 15

cl = makeCluster(no_cluster)
registerDoParallel(cl)

output = foreach( ii=nprev:1, .combine=rbind,
                  .packages = c( "glmnet", "Ohit", "RcppEigen" ) ) %dopar% {
                    sink()
                    rolling.pred(macrodata = dta, nprev = nprev, i = ii,
                                 yar = yar, xar = xar,
                                 eta = 2, q0 = 3, delta0 = 0.5,
                                 training_por = training_por,
                                 grid_size = grid_size,
                                 tune_range = tune_range,
                                 intercept = intercept)
                  }
stopCluster(cl)

Y = dta[,1]
n = length(Y)
acc_sq_err = matrix(0, nrow = nprev, ncol = 7)
acc_abs_err = matrix(0, nrow = nprev, ncol = 7)
error = matrix(0, nrow = nprev, ncol = 7)

for (i in nprev:1) {
  error[(nprev - i + 1),] = Y[n - i + 1] - output[(nprev - i + 1), 1:7]
}
for (t in nprev:1) {
  acc_sq_err[(nprev - t + 1), ] = colSums((error[1:(nprev - t + 1),,drop = F])^2)
  acc_abs_err[(nprev - t + 1), ] = colSums(abs(error[1:(nprev - t + 1),,drop = F]))
}

# Report
cat("Housing starts:", "\n")
cat("y_ar", yar, "x_ar", xar,"\n")
cat("range:", tune_range, "grid_size", grid_size, "prop", training_por, 
    "intercept", intercept, "\n")
cat(c("FHTD", "OGA-3", "AR-OGA", "Lasso", "ALasso", "AR-ALasso", "AR-AIC"), "\n")
# 1. RMSE
round(sqrt( tail( acc_sq_err, 1) / nprev), 4)
round(sqrt(tail( acc_sq_err, 1) / nprev) / sqrt(tail( acc_sq_err, 1) / nprev)[1], 4)

# 2. Median absolute error
round(boxplot(abs(error), plot = F)$stats[3,], 4)
round(boxplot(abs(error), plot = F)$stats[3,] / boxplot(abs(error), plot = F)$stats[3,][1], 4)

# 3. Error analysis
pdf("lines_test3.pdf", width = 12, height = 8)
plot(acc_sq_err[,1], type = "l", ylim = c(0, max(acc_sq_err)), lwd = 1.5, 
     xlab = "rolling windows", ylab = "accumulated squared error", bty = "n")
pal = rainbow(6)
for (i in 2:7) {
  lines(acc_sq_err[,i], col = pal[i - 1])
}
legend("topleft", legend = c("FHTD", "OGA", "arOGA", "Lasso", "ALasso", 
                             "ARALasso", "AR"),
       col = c("Black", pal), lty = rep(1, 7))
dev.off()

# 4. Variable selection analysis
temp = colSums(output)[-c(1:7)]
temp2 = matrix(temp, ncol = 7)

pdf("barplot_test3.pdf", width = 12, height = 6)
barplot(temp2 / nprev, beside = TRUE, ylab = "frequency",
        names.arg = c("FHTD", "OGA-3", "AR-OGA-3", "LASSO", "ALasso",
                      "AR-ALasso", "AR-AIC"))
dev.off()

################################################################################

nprev = 216
intercept = FALSE

################################################################################

no_cluster = 15

cl = makeCluster(no_cluster)
registerDoParallel(cl)

output = foreach( ii=nprev:1, .combine=rbind,
                  .packages = c( "glmnet", "Ohit", "RcppEigen" ) ) %dopar% {
                    sink()
                    rolling.pred(macrodata = dta, nprev = nprev, i = ii,
                                 yar = yar, xar = xar,
                                 eta = 2, q0 = 3, delta0 = 0.5,
                                 training_por = training_por,
                                 grid_size = grid_size,
                                 tune_range = tune_range,
                                 intercept = intercept)
                  }
stopCluster(cl)

Y = dta[,1]
n = length(Y)
acc_sq_err = matrix(0, nrow = nprev, ncol = 7)
acc_abs_err = matrix(0, nrow = nprev, ncol = 7)
error = matrix(0, nrow = nprev, ncol = 7)

for (i in nprev:1) {
  error[(nprev - i + 1),] = Y[n - i + 1] - output[(nprev - i + 1), 1:7]
}
for (t in nprev:1) {
  acc_sq_err[(nprev - t + 1), ] = colSums((error[1:(nprev - t + 1),,drop = F])^2)
  acc_abs_err[(nprev - t + 1), ] = colSums(abs(error[1:(nprev - t + 1),,drop = F]))
}

# Report
cat("Housing starts:", "\n")
cat("y_ar", yar, "x_ar", xar,"\n")
cat("range:", tune_range, "grid_size", grid_size, "prop", training_por, 
    "intercept", intercept, "\n")
cat(c("FHTD", "OGA-3", "AR-OGA", "Lasso", "ALasso", "AR-ALasso", "AR-AIC"), "\n")
# 1. RMSE
round(sqrt( tail( acc_sq_err, 1) / nprev), 4)
round(sqrt(tail( acc_sq_err, 1) / nprev) / sqrt(tail( acc_sq_err, 1) / nprev)[1], 4)

# 2. Median absolute error
round(boxplot(abs(error), plot = F)$stats[3,], 4)
round(boxplot(abs(error), plot = F)$stats[3,] / boxplot(abs(error), plot = F)$stats[3,][1], 4)

# 3. Error analysis
pdf("lines_test4.pdf", width = 12, height = 8)
plot(acc_sq_err[,1], type = "l", ylim = c(0, max(acc_sq_err)), lwd = 1.5, 
     xlab = "rolling windows", ylab = "accumulated squared error", bty = "n")
pal = rainbow(6)
for (i in 2:7) {
  lines(acc_sq_err[,i], col = pal[i - 1])
}
legend("topleft", legend = c("FHTD", "OGA", "arOGA", "Lasso", "ALasso", 
                             "ARALasso", "AR"),
       col = c("Black", pal), lty = rep(1, 7))
dev.off()

# 4. Variable selection analysis
temp = colSums(output)[-c(1:7)]
temp2 = matrix(temp, ncol = 7)

pdf("barplot_test4.pdf", width = 12, height = 6)
barplot(temp2 / nprev, beside = TRUE, ylab = "frequency",
        names.arg = c("FHTD", "OGA-3", "AR-OGA-3", "LASSO", "ALasso",
                      "AR-ALasso", "AR-AIC"))
dev.off()

################################################################################

rm(list=ls())
quit()
