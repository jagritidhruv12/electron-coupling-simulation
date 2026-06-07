library(minpack.lm)
#install.packages("minpack.lm")
library("Hmisc")
#library(Rcpp)


peak_detection <- function(data, nleft, nright, height_threshold){
  index_arr <- -1
  # print("Start peak detection!")
  # print(paste("Length data:", length(data)))
  
  # General strategy: First find all local maxima. 
  # Then check that these are not just local fluctuations but rather actual peaks by
  # requiring that the maximum is at least height_threshold higher than nleft or
  # nright. height_threshold = 1.2 means at least 20% higher
  for(i in (nleft+1):(length(data)-nright)){
    slope_found_left <- TRUE
    for(j in 0:(nleft-1) ){
      # print(paste("left i =",i,"    j =",j))
      # print(paste("Compare ", i-j, "with", i-j-1))
      if(data[i-j] <= data[i-j-1]){
        slope_found_left <- FALSE
      }
    }
    
    slope_found_right <- TRUE
    for(j in 0:(nright-1) ){
      # print(paste("right i =",i,"    j =",j))
      if(data[i+j] <= data[i+j+1] ){
        slope_found_right <- FALSE
      }
    }
    
    if(slope_found_left & slope_found_right){
      if(data[i] > data[i-nleft]*height_threshold &
         data[i] > data[i+nright]*height_threshold){
        index_arr <- c(index_arr,i)
        # print("Peak found!")
      }
    }
  }
  index_arr <- index_arr[-1]
  return(index_arr)
}

find2Peaks <- function(mydata, start_n, additional_error_line, out){
  k <- start_n
  # print(paste("Start find 2 peaks!, k =",k ))
  peaks <- peak_detection(mydata, k, k)
  length_peaks <- length(peaks)
  # print(paste(out, "PeaksA:", length(peaks) ))
  while(length(peaks) != 2 & k>0){
    length_peaks_prev <- length(peaks)
    # is needed to make sure length does not jump between 0/1 and >2 while making only one k-step
    # <<- for global namespace!
    peaks <<- peak_detection(mydata, k,k)
    # print(paste("length peaks =", length(peaks)))
    
    # too many peaks => make condition for peaks more strict 
    if(length(peaks) > 2 & length_peaks_prev != 1){
      k <- k+1 
      # print(paste(out, "k+1 called"))
    }  # too few peaks => make conditions for peaks less strict
    else if(length(peaks) < 2 & length_peaks_prev <= 2){
      k <- k-1    
      # print(paste(out, "k-1 called, k =",k, "length peaks:", length(peaks) ))
    }
    else if(length(peaks) > 2 & length_peaks_prev < 2){
      print(paste(out, "weird case"))
      # this case is if it jumps between e.g. length 1 and 3. Take highest peaks then
      
      index_firstmax <- which( mydata[peaks] == max(mydata[peaks]) )
      print(paste("index firstmax:", index_firstmax))
      
      index_secondmax <- which(mydata[peaks] == 
                                 max( mydata[peaks][mydata[peaks]!=max(mydata[peaks])] ) )
      print(paste("index seoncmdax:", index_secondmax))
      peaks <- peaks[c(index_firstmax, index_secondmax)]
    }
  }
  
  if(k == 0 | length(peaks) != 2){
    if(length(peak_detection(mydata, start_n, start_n)) > length(peaks)){
      peaks <- peak_detection(mydata, start_n, start_n)
    }
    
    print(paste("Could only find", length(peaks),  "out of 2 peaks in FFT:",
                additional_error_line, out) )
    if(length(peaks) == 0){
      peaks <- c(NA,NA)
    }else if(length(peaks) == 1){
      peaks <- c(peaks, NA)
    }
  }
  # print(paste(out, "PeaksB:", length(peaks) ))
  return(peaks)
}




#x_vals_divisor divides the x-values of the fitted data for plotting
#mult_dp is multiplier for plotting lines
fit_Lorentz <- function(tab, span, a_ini, x0_ini, w_ini, y0_ini, plotLinesBool = TRUE,
                        colLine = "black", x_vals_divisor = 1, mult_dp = 200){
  
  
  mypos <- which( (tab[,1]-x0_ini) == min( abs(tab[,1] - x0_ini) ) )
  tab <- tab[c(max(mypos - span/2, 0):min(mypos + span/2, nrow(tab)) ),]
  
  colnames(tab) <- c("x","y")
  
  res <- try(nlsLM( y ~ y0 + (2*a/pi)*w/(4*(x-x0)^2+w^2), 
                    start=c(a=a_ini, x0=x0_ini, w=w_ini, y0=y0_ini),
                    data = tab,  control = list(maxiter = 1024)))
  params <- cbind(coef(res)[[1]], coef(res)[[2]], coef(res)[[3]], 
                  coef(res)[[4]])
  
  eq = function(x, a, x0, w, y0){y0 + (2*a/pi)*w/(4*(x-x0)^2+w^2)  }
  
  
  x_vals_forFit <- seq((floor(min(tab[,1]))*mult_dp),(ceiling(max(tab[,1]))*mult_dp), 
                       length.out=mult_dp*(length(tab[,1])-1)+1)/mult_dp
  fitted_y <- eq(x_vals_forFit,params[1],params[2], params[3], params[4])
  
  if(plotLinesBool == TRUE){
    lines(x_vals_forFit/x_vals_divisor, fitted_y, type='l', col=colLine, lwd=3)
  }
  xc_fit_error <- summary(res)$coefficients[2,2]
  w_fit_error <- summary(res)$coefficients[3,2]
  
  return(c(params, xc_fit_error, w_fit_error))
  
}

#x_vals_divisor divides the x-values of the fitted data for plotting
#mult_dp is multiplier for plotting lines
fit_Lorentz2 <- function(tab, span, a_ini, x0_ini, b_ini, y0_ini = 0, plotLinesBool = TRUE,
                         colLine = "red", x_vals_divisor = 1, mult_dp = 200){
  
  colnames(tab) <- c("x","y")
  
  # Take index of the frequency closest to x0_ini and cut on span data points
  f_lower <- x0_ini - span
  f_higher <- x0_ini + span
  tab <- tab[which(tab$x > f_lower & tab$x < f_higher ),]
  
  res <- try(nlsLM( y ~ a/((x-x0)^2+b^2) + y0,
                    start=c(a=a_ini, x0=x0_ini, b=b_ini, y0=y0_ini),
                    lower=c(0, 0, 0, -Inf),
                    data = tab,  control = list(maxiter = 1024)))
  params <- cbind(coef(res)[[1]], coef(res)[[2]], coef(res)[[3]], coef(res)[[4]])
  
  eq = function(x, a, x0, b, y0){ a/((x-x0)^2+b^2) + y0 }
  
  
  x_vals_forFit <- seq((floor(min(tab[,1]))*mult_dp),(ceiling(max(tab[,1]))*mult_dp), 
                       length.out=mult_dp*(length(tab[,1])-1)+1)/mult_dp
  fitted_y <- eq(x_vals_forFit,params[1],params[2], params[3], params[4])
  print(x_vals_forFit[which(fitted_y == max(fitted_y)) ])
  
  if(plotLinesBool == TRUE){
    lines(x_vals_forFit/x_vals_divisor, fitted_y, type='l', col=colLine, lwd=3)
  }
  xc_fit_error <- summary(res)$coefficients[2,2]
  b_fit_error <- summary(res)$coefficients[3,2]
  
  return(c(params, xc_fit_error, b_fit_error))
  
}


fit_Lorentz3 <- function(tab, span, a_ini, x0_ini, b_ini, y0_ini = 0, plotLinesBool = TRUE,
                         colLine = "red", x_vals_divisor = 1, mult_dp = 200){
  
  colnames(tab) <- c("x","y")
  
  # Take index of the frequency closest to x0_ini and cut on span data points
  f_lower <- x0_ini - span
  f_higher <- x0_ini + span
  tab <- tab[which(tab$x > f_lower & tab$x < f_higher ),]
  
  res <- try(nlsLM( y ~ a/(4*x0^2) * 1/((x-x0)^2 + b^2/4) + y0,
                    start=c(a=a_ini, x0=x0_ini, b=b_ini, y0=y0_ini),
                    lower=c(0, 0, 0, -Inf),
                    data = tab,  control = list(maxiter = 1024)))
  params <- cbind(coef(res)[[1]], coef(res)[[2]], coef(res)[[3]], coef(res)[[4]])
  
  eq = function(x, a, x0, b, y0){ a/(4*x0^2) * 1/((x-x0)^2 + b^2/4) + y0 }
  
  
  x_vals_forFit <- seq((floor(min(tab[,1]))*mult_dp),(ceiling(max(tab[,1]))*mult_dp), 
                       length.out=mult_dp*(length(tab[,1])-1)+1)/mult_dp
  fitted_y <- eq(x_vals_forFit,params[1],params[2], params[3], params[4])
  print(x_vals_forFit[which(fitted_y == max(fitted_y)) ])
  
  if(plotLinesBool == TRUE){
    lines(x_vals_forFit/x_vals_divisor, fitted_y, type='l', col=colLine, lwd=3)
  }
  xc_fit_error <- summary(res)$coefficients[2,2]
  b_fit_error <- summary(res)$coefficients[3,2]
  
  return(c(params, xc_fit_error, b_fit_error))
  
}

eqLorentzSqrt <- function(x, a, x0, b, y0){ sqrt(a/((x-x0)^2+b^2)) + y0 }

fitLorentzSqrt <- function(tab, span, a_ini, x0_ini, b_ini, y0_ini, plotLinesBool = TRUE,
                           colLine = "black", x_vals_divisor = 1, mult_dp = 200){
  
  # Take index of the frequency closest to x0_ini and cut on span data points
  mypos <- which( (tab[,1]-x0_ini) == min( abs(tab[,1] - x0_ini) ) )
  tab <- tab[c(max(mypos - span/2, 0):min(mypos + span/2, nrow(tab)) ),]
  
  colnames(tab) <- c("x","y")
  
  res <- try(nlsLM( y ~ sqrt(a/((x-x0)^2+b^2)) + y0,
                    start=c(a=a_ini, x0=x0_ini, b=b_ini, y0=y0_ini),
                    lower=c(0, 0, 0, -Inf),
                    data = tab,  control = list(maxiter = 1024)))
  params <- cbind(coef(res)[[1]], coef(res)[[2]], coef(res)[[3]], coef(res)[[4]])
  
  eq = function(x, a, x0, b, y0){ sqrt(a/((x-x0)^2+b^2)) + y0 }
  
  
  x_vals_forFit <- seq((floor(min(tab[,1]))*mult_dp),(ceiling(max(tab[,1]))*mult_dp), 
                       length.out=mult_dp*(length(tab[,1])-1)+1)/mult_dp
  fitted_y <- eq(x_vals_forFit,params[1],params[2], params[3], params[4])
  print(x_vals_forFit[which(fitted_y == max(fitted_y)) ])
  
  if(plotLinesBool == TRUE){
    lines(x_vals_forFit/x_vals_divisor, fitted_y, type='l', col=colLine, lwd=3)
  }
  xc_fit_error <- summary(res)$coefficients[2,2]
  b_fit_error <- summary(res)$coefficients[3,2]
  
  return(c(params, xc_fit_error, b_fit_error))
  
}


fitSingleDip <- function(tab, N_particles, span_fit, plotLinesBool = TRUE){
  
  f_lower_fit <- 626951 - span_fit/2
  f_higher_fit <- 626951 + span_fit/2
  tab2 <- tab[which(tab[,1] > f_lower_fit & tab[,1] < f_higher_fit ),]
  if(nrow(tab2) < 10){
    stop("Data has less than 10 rows left :(")
  }
  colnames(tab2) <- c("f", "y")
  
  n0start <- max(tab2[,2])
  n1start <- -3
  tab3 <- tab[which(tab[,1] > 626951-10 & tab[,1] < 626951+10 ),]
  f0start <- tab3[which(tab3[,2] == min(tab3[,2]))[1],1]
  f1start <- f0start
  df0start <- 200
  df1start <- 0.1*N_particles
  # https://www.wolframalpha.com/input/?i=is+4%28%28f-f1%29%2Fd1%29%5E2+%2F%284%28%28f-f1%29%2Fd1%29%5E2+%2B+%284*%28f-f0%29%2Fd0+*+%28f-f1%29%2Fd1+-+1%29%5E2%29+%3D%3D+4%28f-f1%29%5E2%2F%284%28f-f1%29%5E2+%2B+%284%28f-f0%29%2Fd0+*+%28f-f1%29+-+d1%29%5E2%29
  res <- try(nlsLM( y ~ n0+ n1*4*(f-f1)^2/(4*(f-f1)^2+(4*((f-f0)/df0)*(f-f1)-df1)^2), 
                    start=c(n0=n0start, n1=n1start, f0=f0start, df0=df0start, f1=f1start, df1=df1start),
                    data = tab2, lower=c(-Inf,-Inf,0,0,0,0), control = list(maxiter = 1024)))
  
  params <- cbind(coef(res)[[1]],coef(res)[[2]],coef(res)[[3]],coef(res)[[4]],coef(res)[[5]],coef(res)[[6]])
  
  eq = function(f,n0,n1,f0,df0,f1,df1){n0+n1*4*(f-f1)^2/(4*(f-f1)^2+(4*(f-f0)/df0*(f-f1)-df1)^2)}
  if(plotLinesBool){
    lines(tab2[,1], eq(tab2[,1], params[1], params[2], params[3], params[4], params[5], params[6]), 
          type='l', col="red", lwd=2)
  }
  return(params)
}

eqSingleDipSqrt <- function(f,n0,n1,f0,df0,f1,df1){
  n0 + sqrt(n1*4*(f-f1)^2/(4*(f-f1)^2+(4*(f-f0)/df0*(f-f1)-df1)^2))
}


fitSingleDipSqrt <- function(tab, span_fit, n0ini = 0, n1ini, f0ini = 626951.5 ,df0ini = 20, f1ini = 626951.5, df1ini = 5,
                             plotLinesBool = TRUE, plotStartingValues = FALSE){
  
  f_lower_fit <- 626951 - span_fit/2
  f_higher_fit <- 626951 + span_fit/2
  tab2 <- tab[which(tab[,1] > f_lower_fit & tab[,1] < f_higher_fit ),]
  if(nrow(tab2) < 10){
    stop("Data has less than 10 rows left :(")
  }
  colnames(tab2) <- c("f", "y")
  
  # n0start <- 0
  # n1start <- max(tab2$y)^2
  # tab3 <- tab[which(tab[,1] > 626951-10 & tab[,1] < 626951+10 ),]
  # f0start <- tab3[which(tab3[,2] == min(tab3[,2]))[1],1]
  # f1start <- f0start
  # df0start <- 20
  # df1start <- 5
  res <- try(nlsLM( y ~ n0+ sqrt(n1*4*(f-f1)^2/(4*(f-f1)^2+(4*((f-f0)/df0)*(f-f1)-df1)^2)), 
                    start=c(n0=n0ini, n1=n1ini, f0=f0ini, df0=df0ini, f1=f1ini, df1=df1ini),
                    data = tab2, lower=c(-Inf,-Inf,0,0,0,0), control = list(maxiter = 1024)))
  
  params <- cbind(coef(res)[[1]],coef(res)[[2]],coef(res)[[3]],coef(res)[[4]],coef(res)[[5]],coef(res)[[6]])
  
  eq = function(f,n0,n1,f0,df0,f1,df1){n0 + sqrt(n1*4*(f-f1)^2/(4*(f-f1)^2+(4*(f-f0)/df0*(f-f1)-df1)^2)) }
  
  if(plotStartingValues){
    lines(tab2$f, eq(tab2$f, n0ini, n1ini, f0ini, df0ini, f1ini, df1ini), col="blue")
  }
  
  if(plotLinesBool){
    lines(tab2[,1], eq(tab2[,1], params[1], params[2], params[3], params[4], params[5], params[6]), 
          type='l', col="red", lwd=2)
  }
  return(params)
}


fitSingleDipSqrtlog10 <- function(tab, span_fit, n0ini = 0, n1ini, f0ini = 626951.5 ,df0ini = 20, f1ini = 626951.5, df1ini = 5,
                             plotLinesBool = TRUE, plotStartingValues = FALSE){
  
  f_lower_fit <- 626951 - span_fit/2
  f_higher_fit <- 626951 + span_fit/2
  tab2 <- tab[which(tab[,1] > f_lower_fit & tab[,1] < f_higher_fit ),]
  if(nrow(tab2) < 10){
    stop("Data has less than 10 rows left :(")
  }
  colnames(tab2) <- c("f", "y")
  
  # n0start <- 0
  # n1start <- max(tab2$y)^2
  # tab3 <- tab[which(tab[,1] > 626951-10 & tab[,1] < 626951+10 ),]
  # f0start <- tab3[which(tab3[,2] == min(tab3[,2]))[1],1]
  # f1start <- f0start
  # df0start <- 20
  # df1start <- 5
  res <- try(nlsLM( y ~ n0+ sqrt(n1*4*(f-f1)^2/(4*(f-f1)^2+(4*((f-f0)/df0)*(f-f1)-df1)^2)) , 
                    start=c(n0=n0ini, n1=n1ini, f0=f0ini, df0=df0ini, f1=f1ini, df1=df1ini),
                    data = tab2, lower=c(-Inf,-Inf,0,0,0,0), control = list(maxiter = 1024)))
  
  params <- cbind(coef(res)[[1]],coef(res)[[2]],coef(res)[[3]],coef(res)[[4]],coef(res)[[5]],coef(res)[[6]])
  
  eq = function(f,n0,n1,f0,df0,f1,df1){n0 + sqrt(n1*4*(f-f1)^2/(4*(f-f1)^2+(4*(f-f0)/df0*(f-f1)-df1)^2)) }
  
  if(plotStartingValues){
    lines(tab2$f, eq(tab2$f, n0start, n1start, f0start, df0start, f1start, df1start), col="blue")
  }
  
  if(plotLinesBool){
    lines(tab2[,1], eq(tab2[,1], params[1], params[2], params[3], params[4], params[5], params[6]), 
          type='l', col="red", lwd=2)
  }
  return(params)
}

fitDoubleDip <- function(tab, span_fit,  f1start, f2start, plotLinesBool = TRUE,
                         N_particles1 = 1, N_particles2 = 1){
  
  result = tryCatch({
    f_lower_fit <- 626951 - span_fit/2
    f_higher_fit <- 626951 + span_fit/2
    tab2 <- tab[which(tab[,1] > f_lower_fit & tab[,1] < f_higher_fit ),]
    if(nrow(tab2) < 10){
      stop("Data has less than 10 rows left :(")
    }
    colnames(tab2) <- c("f", "y")
    
    n0start <- max(tab2[,2])
    n1start <- -3
    # tab3 <- tab[which(tab[,1] > 626951-10 & tab[,1] < 626951+10 ),]
    f0start <- 626951
    df0start <- 200
    df1start <- 0.6*N_particles1
    df2start <- 0.6*N_particles2
    
    res <- nlsLM( y ~ n0+n1*4*((f-f1)/df1)^2*((f-f2)/df2)^2/(4*((f-f1)/df1)^2*((f-f2)/df2)^2+(4*((f-f0)/df0)*((f-f1)/df1) * ((f-f2)/df2) - (f-f1)/df1 - (f-f2)/df2  )^2), 
                  start=c(n0=n0start, n1=n1start, f0=f0start, df0=df0start, f1=f1start, df1=df1start, f2=f2start, df2 = df2start),
                  data = tab2,  control = list(maxiter = 1024))
    
    params <- cbind(coef(res)[[1]],coef(res)[[2]],coef(res)[[3]],coef(res)[[4]],coef(res)[[5]],coef(res)[[6]], coef(res)[[7]], coef(res)[[8]])
    
    eq = function(f,n0,n1,f0,df0,f1,df1,f2,df2){n0+n1*4*((f-f1)/df1)^2*((f-f2)/df2)^2/(4*((f-f1)/df1)^2*((f-f2)/df2)^2+(4*((f-f0)/df0)*((f-f1)/df1) * ((f-f2)/df2) - (f-f1)/df1 - (f-f2)/df2  )^2)}
    if(plotLinesBool){
      lines(tab2[,1], eq(tab2[,1], params[1], params[2], params[3], params[4], params[5], params[6], params[7], params[8]), 
            type='l', col="red", lwd=1)
    }
  }, warning = function(warning_condition) {
    # warning-handler-code
  }, error = function(error_condition) {
    params <- rep(-1, times=8)
  }, finally={
    return(params)
  })
  
}

