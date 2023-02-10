#Utility functions for mhmm

#Outputs the map for chosen covariates. 
# Beta: the full covariance structure
# sel_cov: vector of TRUE/FALSE. TRUE = include covariate in 
# row corresponding to position in sel_cov
Beta_map <- function(Beta, sel_cov){
  
  map.matrix <- Beta
  map.matrix[sel_cov, ] <- seq(sum(sel_cov)*dim(Beta)[2])
  map.matrix[!sel_cov, ] <- NA
  return(factor(map.matrix))
  
}

# Function to transform natural parameters to working

delta.n2w <- function(m, delta){
  
  foo <- log(delta/delta[1])
  tdelta <- as.vector(tail(foo, m - 1))
  return(tdelta) 
}

delta.w2n <- function(m, tdelta){
  
  # set first element to one and fill in the last m - 1 elements with working parameters and take exp
  delta <- c(1, exp(tdelta))
  
  # normalize
  delta = delta/sum(delta)
  
  return(delta)
}
#Function to transform natural parameters to working (C ++ code as basis then add 1 to elements)
gamma.n2w <- function(m, gamma){
  
  foo <- log(gamma/diag(gamma))
  
  tgamma <- seq(m*(m-1))
  
  for(idx in 0:((m*(m - 1)/2 - 1))){
    i <-  1 + as.integer((- 1 + sqrt(8*idx + 1))/2)
    j <- idx - (i - 1)*i/2
    #Adjust c++ code
    i <- i + 1
    j <- j + 1
    idx <- idx + 1
    
    tgamma[idx] <- foo[i,j]
    tgamma[m*(m-1)/2 + idx] <- foo[j,i]
    
  }
  
  return(tgamma)
}

#Function tranforming working parameters to natural parameters (C ++ code as basis then add 1 to elements)
gamma.w2n <- function(m, tgamma){
  
  # Construct m x m identity matrix
  gamma <- diag(m)
  
  #Fill offdiagonal elements with working parameters colwise:  
  for(idx in 0:((m*(m - 1)/2 - 1))){
    
    i <-  1 + as.integer((- 1 + sqrt(8*idx + 1))/2)
    j <- idx - (i - 1)*i/2
    
    #Adjust c++ code
    i <- i + 1
    j <- j + 1
    idx <- idx + 1
    
    #fill gamma according to mapping and take exp()  
    gamma[i,j] = exp(tgamma[idx])
    gamma[j,i] = exp(tgamma[m*(m-1)/2 + idx])
    
  }
  
  #Normalize each row
  gamma <- gamma/apply(gamma, 1, sum)
  
  return(gamma);
}


# Computes smoothing probabilities and decode
hmm_decode <- function(MakeADFun.obj)
{
  #Get objects we need from MakeADFun object
  adreps <- MakeADFun.obj$report(MakeADFun.obj$env$last.par.best) #Retrieve the objects at ML value
  lambda_mat <- adreps$lambda_mat
  filt_prob <- adreps$filt_prob
  pred_prob <- adreps$pred_prob
  gamma_lag <- adreps$gamma_lag
  delta_lag <- adreps$delta_lag
  m <- adreps$m
  T <- adreps$T
  mstar <- dim(gamma_lag)[1]
  perm <- adreps$perm
  eta0 <- adreps$eta0
  y0 <- adreps$y0
  y <- adreps$y
  a <- adreps$a
  b <- adreps$b
  d <- adreps$d
  
  # Compute smoothing probabilities for S^*t
  smooth_prob <- matrix(0, nrow = T, ncol = mstar)
  smooth_prob[T, ] <- filt_prob[T, ]
  for(t in (T - 1):1){
    for(i in 1:mstar){
      smooth_prob[t, i] <- filt_prob[t, i]*sum(smooth_prob[t + 1, ]*gamma_lag[i, ]/(pred_prob[t, ] + 1e-10))   
    }
  }
  
  # Compute corresponding smoothing probabilities of S^t
  stateprobs <- matrix(NA, nrow = T, ncol = m)
  selvec <- perm[, 2] # Define which vector
  for(i in 1:m){
    stateprobs[ ,i] <- as.vector(rowSums(smooth_prob[ ,selvec == (i - 1)])) # i - 1 here because of compatibility with perm
  }
  
  # Local decoding using smoothing probs
  ldecode <- rep(NA, T)
  for (i in 1:T) ldecode[i] <- which.max(stateprobs[i, ])
  
  
  
  ### In sample prediction - standard way ###
  
  #Define P(S^*_t | Omega_{t - 1}) matrix
  stateprobs_pred_lag <- matrix(NA, nrow = T, ncol = m^2)
  stateprobs_pred_lag[1, ] <- delta_lag                          # Let P(S_1 | Omega_0) = stationary distribution
  stateprobs_pred_lag[2:T, ] <- pred_prob[1:(T - 1), ]           # The algorithm produces P(S^*_(t + 1) | Omega_{t}) at time t
  
  #Find corresponding values of P(S^_t | Omega_{t - 1})
  stateprobs_pred <- matrix(NA, nrow = T, ncol = m)
  for(i in 1:m){
    stateprobs_pred[ ,i] <- as.vector(rowSums(stateprobs_pred_lag[ ,selvec == (i - 1)])) # i - 1 here because of compatibility with perm
  }
  
  # Local decoding 
  ldecode_pred <- rep(NA, T)
  for (i in 1:T) ldecode_pred[i] <- which.max(stateprobs_pred[i, ])
  
  #Predicted intensity
  lambda_pred <- rowSums(stateprobs_pred_lag*lambda_mat)
  
  # Predicted intensity - using smoothing probabilities ##
  lambda_pred_smooth <- rowSums(smooth_prob*lambda_mat)
  
  # Output as list
  list(stateprobs = stateprobs, ldecode = ldecode, lambda_pred_smooth = lambda_pred_smooth, stateprobs_pred = stateprobs_pred, ldecode_pred = ldecode_pred, lambda_pred = lambda_pred)
}



# Nice publication-worthy theme
theme_Publication <- function(base_size=12) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size = base_size)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(size = rel(0.8)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "top",
            legend.direction = "horizontal",
            legend.key.size= unit(0.5, "cm"),
            legend.spacing = unit(0, "cm"),
            legend.title = element_text(),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text()
    ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}


# QQplot with lines
gg_qq <- function(x, distribution = "norm", ..., line.estimate = NULL, conf = 0.95,
                  labels = names(x)){
  q.function <- eval(parse(text = paste0("q", distribution)))
  d.function <- eval(parse(text = paste0("d", distribution)))
  x <- na.omit(x)
  ord <- order(x)
  n <- length(x)
  P <- ppoints(length(x))
  df <- data.frame(ord.x = x[ord], z = q.function(P, ...))
  
  if(is.null(line.estimate)){
    Q.x <- quantile(df$ord.x, c(0.25, 0.75))
    Q.z <- q.function(c(0.25, 0.75), ...)
    b <- diff(Q.x)/diff(Q.z)
    coef <- c(Q.x[1] - b * Q.z[1], b)
  } else {
    coef <- coef(line.estimate(ord.x ~ z))
  }
  
  zz <- qnorm(1 - (1 - conf)/2)
  SE <- (coef[2]/d.function(df$z)) * sqrt(P * (1 - P)/n)
  fit.value <- coef[1] + coef[2] * df$z
  df$upper <- fit.value + zz * SE
  df$lower <- fit.value - zz * SE
  
  if(!is.null(labels)){ 
    df$label <- ifelse(df$ord.x > df$upper | df$ord.x < df$lower, labels[ord],"")
  }
  
  p <- ggplot(df, aes(x=z, y=ord.x)) +
    geom_point() + 
    geom_abline(intercept = coef[1], slope = coef[2]) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha=0.2) 
  if(!is.null(labels)) p <- p + geom_text( aes(label = label))
  print(p)
  coef
}


### Mean checker! ###########
mean.check <- function (a, b, d, gamma, n = 5000){
  
  delta <- solve(t(diag(2) - gamma + 1), rep (1,2)) #assume stationary
  eta0 <- delta[1]*d[1]/(1 - a[1] - b[1]) + delta[2]*d[2]/(1 - a[2] - b[2])
  y0 <- exp(eta0)
  
  # Seed
  set.seed(1)
  # set.seed(9) #OK
  
  ##Sample time series
  eta <- eta1 <- eta2 <-  y <- y1 <- y2 <- rep(NA,n)
  #sample init state from delta
  state0 = sample(c(1,2), size = 1, prob = delta)
  #sample chain
  outs <- rmarkovchain(n = n, object = new("markovchain", states = c("1", "2"),  
                                           transitionMatrix = gamma), what = "list", t0 = state0)
  # t=1
  if(state0 == "1"){eta[1] <- d[1] + a[1]*eta0 + b[1]*log(y0 + 1)}
  if(state0 == "2"){eta[1] <- d[2] + a[2]*eta0 + b[2]*log(y0 + 1)}
  
  y[1] <- rpois(1, lambda = exp(eta[1]))
  
  #t=2,...,n
  for(i in 2:n){
    if(outs[i] == "1"){eta[i] <- d[1] + a[1]*eta[i - 1] + b[1]*log(y[i - 1] + 1)}
    if(outs[i] == "2"){eta[i] <- d[2] + a[2]*eta[i - 1] + b[2]*log(y[i - 1] + 1)}
    y[i] <- rpois(1, lambda = exp(eta[i]))
  }
  
  mean.state1 <- mean(y[outs == "1"])
  mean.state2 <- mean(y[outs == "2"])
  
  means <- c(mean.state1, mean.state2)
  names(means) <- c("mean state 1", "mean state 2")
  
  # Marginal
  # t=1
  eta1[1] <- d[1] + a[1]*eta0 + b[1]*log(y0 + 1)
  eta2[1] <- d[2] + a[2]*eta0 + b[2]*log(y0 + 1)
  
  y1[1] <- rpois(1, lambda = exp(eta1[1]))
  y2[1] <- rpois(1, lambda = exp(eta2[1]))
  
  for(i in 2:n){
    eta1[i] <- d[1] + a[1]*eta1[i - 1] + b[1]*log(y1[i - 1] + 1)
    eta2[i] <- d[2] + a[2]*eta2[i - 1] + b[2]*log(y2[i - 1] + 1)
    y1[i] <- rpois(1, lambda = exp(eta1[i]))
    y2[i] <- rpois(1, lambda = exp(eta2[i]))
  }
  
  means.marginal <- c(mean(y1), mean(y2))
  names.means <- c("marginal mean state 1", "marginal mean state 2")
  return(rbind(means, means.marginal))
  
}


# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


#######  function generating linear constraints for m = 2,3  #########################################
# used with constrOtim (stats)
constr.fun <- function(m, ncovpar = 2, eps = 1e-7){
    if(m == 2){  
      # Linear constraints when m = 2, in the formula ui %*% theta - ci >= 0, ci = (cci + eps)
      ui <- matrix(0, nrow = 10, ncol = 8 +  + ncovpar)
      ui[1, ] <- c(-1, rep(0, 7 + ncovpar))           # for a1 < 1
      ui[2, ] <- c(1, rep(0, 7 + ncovpar))            # for a1 > 0
      ui[3, ] <- c(0, -1, rep(0, 6 + ncovpar))        # for a2 < 1 
      ui[4, ] <- c(0, 1, rep(0, 6 + ncovpar))         # for a2 > 0 
      ui[5, ] <- c(0, 0, -1, rep(0, 5 + ncovpar))     # for b1 < 1
      ui[6, ] <- c(0, 0, 1, rep(0, 5 + ncovpar))      # for b1 > 0
      ui[7, ] <- c(0, 0, 0, -1, rep(0, 4 + ncovpar))  # for b2 < 1
      ui[8, ] <- c(0, 0, 0, 1, rep(0, 4 + ncovpar))   # for b2 > 0
      ui[9, ] <- c(-1, 0, -1, rep(0, 5 + ncovpar))    # for a1 + b1 < 1
      ui[10, ] <- c(0, -1, 0, -1, rep(0,4 + ncovpar)) # for a2 + b2 < 1
      
      cci <- c(-1, 0, -1, 0, -1, 0, -1, 0, -1, -1) 
      ci <- cci + eps
    }# end of if m = 2
    if(m == 3){  
      # Linear constraints when m = 2, in the formula ui %*% theta - ci >= 0, ci = (cci + eps)
      ui <- matrix(0, nrow = 15, ncol = 15)
      ui[1, ] <- c(-1, rep(0, 14 ))          # for a1 < 1
      ui[2, ] <- c(1, rep(0, 14))            # for a1 > 0
      ui[3, ] <- c(0, -1, rep(0, 13))        # for a2 < 1 
      ui[4, ] <- c(0, 1, rep(0, 13))         # for a2 > 0 
      ui[5, ] <- c(0, 0, -1, rep(0, 12))     # for a3 < 1
      ui[6, ] <- c(0, 0, 1, rep(0, 12))      # for a3 > 0
      ui[7, ] <- c(0, 0, 0, -1, rep(0, 11))           # for b1 < 1
      ui[8, ] <- c(0, 0, 0, 1, rep(0, 11))            # for b1 > 0
      ui[9, ] <- c(0, 0, 0, 0, -1, rep(0, 10))        # for b2 < 1
      ui[10, ] <- c(0, 0, 0, 0, 1, rep(0, 10))         # for b2 > 0
      ui[11, ] <- c(0, 0, 0, 0, 0, -1, rep(0, 9))     # for b3 < 1
      ui[12, ] <- c(0, 0, 0, 0, 0, 1, rep(0, 9))      # for b3 > 0
      ui[13, ] <- c(-1, 0, 0, -1, rep(0, 11))         # for a1 + b1 < 1
      ui[14, ] <- c(0, -1, 0, 0, -1, rep(0, 10))       # for a2 + b2 < 1
      ui[15, ] <- c(0, 0, -1, 0, 0, -1, rep(0, 9))    # for a3 + b3 < 1
      
      cci <- c(rep(c(-1, 0), 6), # for the individual 0 < ai,bi < 1
               -1, -1, -1)        # for the ai + bi < 1
      ci <- cci + eps
    } # end of if m = 3
    
    return(list(ui = ui, ci = ci))
  }  
  

# For screening outliers
remove_outliers <- function(x, na.rm = TRUE, probs=c(.25, .75), ...) {
  qnt <- quantile(x, probs=probs, na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
  # y <- y[complete.cases(y)]
}

# For adding labels
add_label <- function(xfrac, yfrac, label, pos = 4, ...) {
  u <- par("usr")
  x <- u[1] + xfrac * (u[2] - u[1])
  y <- u[4] - yfrac * (u[4] - u[3])
  text(x, y, label, pos = pos, ...)
}