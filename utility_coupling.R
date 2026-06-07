library(minpack.lm)
#install.packages("minpack.lm")
library("Hmisc")
library(Rcpp)
source("./utility_fitting.R")

# data <- ps_save$vx
# sample_freq <- 1/(t_step*save_every_xth)
# plottitle <- plottitle

qC2U0_p <- 1.29776307354e-14
mp <- 1.67262192369e-27
w_z <- sqrt(2*qC2U0_p/mp)

bin_mean <- function (vec, every, na.rm = FALSE) {
  n <- length(vec)
  x <- .colMeans(vec, every, n %/% every, na.rm)
  r <- n %% every
  if (r) x <- c(x, mean.default(vec[(n - r + 1):n], na.rm = na.rm))
  x
}


write.table_with_header <- function(x, file, header, ...){
  cat(header, '\n',  file = file)
  write.table(x, file, append = T, ...)
}




setupParallelComputing <- function(number_of_steps, save_every_xth,
                                   paramToScan, number_of_cores){
  
  if(Sys.info()["nodename"] == "QUANTUM19-9" & number_of_cores < 15){
    # stop("Auf Quantum19-9 weniger als 15 cores angestellt! ")
  }
  
  #setup parallel backend to use many processors
  maxcores <- detectCores()
  cores <- 0
  if(number_of_cores < 0 | number_of_cores > maxcores ){
    cores <- maxcores
  }else{
    cores <- number_of_cores
  }
  # cores=detectCores()
  close( file( "./monitorfile.txt", open="w" ) )
  cl <- makeCluster(min(cores, length(paramToScan)), outfile="./monitorfile.txt") #to not overload your computer
  registerDoParallel(cl)
  
  print(paste(min(cores, length(paramToScan)), " core started") )
  # save_every_xth <- 750
  
  # number of columns for your resulting data. 5 for t, zp, vzp, zBe, vzBe
  columnsResultMatrix <- 3
  totalDatapoints <- number_of_steps*columnsResultMatrix/save_every_xth *
    min((cores[1]-1), length(paramToScan))
  # RAM_needed <- 0.4*16/7e8 * totalDatapoints 
  RAM_needed <- totalDatapoints/(15e9*15/750 *5) * 0.45*32
  
  # array mit 15e9 * 15/750 rows und 5 columns entspricht 45% RAM
  
  free_RAM_GB <- as.numeric(gsub("\r","",gsub("FreePhysicalMemory=","",
                                              system('wmic OS get FreePhysicalMemory /Value',
                                                     intern=TRUE)[3])))/1024/1024
  if(free_RAM_GB < RAM_needed ){
    stop(paste("You want to use ", RAM_needed, "GB of your free ", free_RAM_GB, "GB",
               sep=""))
    stopCluster(cl)
  }
  if(free_RAM_GB < RAM_needed*1.2 ){
    message(paste("You will use an estimated ",
                  RAM_needed*100/free_RAM_GB, "% of your free RAM.
                  Do you want to continue? (type y for yes, anything for no)", sep="" ))
    answer <- readline()
    if(answer!="y"){
      stopCluster(cl)
      stop("Execution successfully aborted.")
    }
  }
  # print("Careful, RAM-check disabled!")
  

  clusterExport(cl, list("nlsLM", "rollmean", "SimulationRun", "ceil") )
  return(cl)
  
}


get_spectrum_no_str <- function(i){
  if(i<10){
    spectrum_no_str <- paste("00", toString(i), sep="")
  }else if(i<100){
    spectrum_no_str <- paste("0", toString(i), sep="")
  }else{
    spectrum_no_str <- toString(i)
  }
  spectrum_no_str
}


doSomeChecks <- function(t_step, number_of_steps, save_every_xth, bool_saveCyclotronFFT, fplus_Be){
  sample_freq <- 1/(t_step*save_every_xth)
  if(sample_freq < 650e3){
    stop(paste("Your maximum FFT frequency is", sample_freq, 
               "\nDecrease t_step or save_every_xth to continue!"))
  }
  
  # number_of_steps <- 0.01e9
  # save_every_xth <- 7.5e3
  # check if number of entries is a floating point (has caused errors in the past)
  temp <- number_of_steps/save_every_xth
  if(temp%%1 != 0 ){
    stop("Your number of entries is not an integer! Adjust number_of_steps or save_every_xth!")
  }
  
  if(bool_saveCyclotronFFT & 1/(dt*save_every_xth) < fplus_Be){
    stop(paste("Cyclotron data save flag set to true, but max FFT frequency at ", 1/(dt*save_every_xth) * 1e-6, "MHz"), sep="")
  }
  
}

doSomeChecks_SimPC <- function(t_step, number_of_steps, save_every_xth){
  sample_freq <- 1/(t_step*save_every_xth)
  if(sample_freq < 28e6){
    stop(paste("Your maximum FFT frequency is", sample_freq, 
               "\nDecrease t_step or save_every_xth to continue!"))
  }
  
  # number_of_steps <- 0.01e9
  # save_every_xth <- 7.5e3
  # check if number of entries is a floating point (has caused errors in the past)
  temp <- number_of_steps/save_every_xth
  if(temp%%1 != 0 ){
    stop("Your number of entries is not an integer! Adjust number_of_steps or save_every_xth!")
  }

  # 
  # if(min(length(paramToScan), 15) *number_of_steps/save_every_xth * 1/(15e9*15/750) < 80 )
  
}


getRadialAmplitudes <- function(data, fcenter, number_of_periods){
  
  n <- round(number_of_periods*1/(fcenter*dt))
  n <- 2^ceil(log2(n))
  n <- min(n, nrow(data))
  j <- 1
  amplitudes <- data.frame(t = c(1:floor(nrow(data)/(j*n))), rho = c(1:floor(nrow(data)/(j*n))) )
  while(j*n < nrow(data)){
    tab <- data[c(1:n) + (j-1)*n,]
    tmean <- mean(tab$t)
    tab$f <- seq(0, nrow(tab)-1,by=1)*(1/(dt*save_every_xth*nrow(tab)  ))
    tab$fft <- 2*abs(fft(tab$xBe))/nrow(tab)
    tab <- tab[which(tab$f < (fcenter + 10e3) & tab$f > (fcenter - 10e3)) ,  ]
    # plot(tab$f, tab$fft)
    
    amplitudes[j,] <-c(tmean, sqrt(sum(tab$fft^2)))
    
    j <- j+1
  }
  
  amplitudes
  
}

# lot_and_save_FFT(ps_save$zp + ps_save$zBe, dt, save_every_xth, fz_p-5e3, 20e3, 10,
                 # plotfilename, spectrum_no_str, plottext, filename_data, bool_vert_line = FALSE)
# data <- ps_save$zp + ps_save$zBe
# dt <- dt
# f_center <- fz_p - 5e3
# f_range <- 20e3
# data_points_to_av <- 10
# plotfilename <- plotfilename
# plottext <- ""
# # filename
plot_and_save_FFT <- function(data, dt, save_every_xth, f_center, f_range, data_points_to_av,
                              plotfilename, spectrum_no_str, plottext, filename_data, bool_vert_line = TRUE){
  sample_freq <- 1/(dt*save_every_xth)
  
  
  print("Start calculating fft...")
  # multiplied by 2 because of alias symmetry (see factor of 2 demonstration R code). psd and not the root bc power must be conserved
  # /length(data) bc of energy conservation theorem
  # 1/df to get psd, and df=1/dt is the bandwidth every point "got filled with"
  psd <- 2*abs(fft(data))^2/length(data)* 1/sample_freq
  freq_axis <- seq(0, length(data)-1,by=1)*(sample_freq/length(data))
  
  f_lower <- f_center - f_range
  f_higher <- f_center + f_range
  psd <- psd[which(freq_axis > f_lower & freq_axis < f_higher )]
  freq_axis <- freq_axis[which(freq_axis > f_lower & freq_axis < f_higher )]
  
  # data_points_to_av <- 10
  # psd_U_log <- log10(psd_U)
  if(min(psd) > 0){
    freq_axis_smooth <- rollmean(freq_axis, data_points_to_av)
    psd_log_smooth <- rollmean(log10(psd), data_points_to_av)
    png(plotfilename, width=1080*2, height=720*2, pointsize=20*2)
    plot(freq_axis_smooth, psd_log_smooth, type='l',
         xlab = "Frequency [Hz]", ylab= "log10(data) [V]")
    mtext(plottext )
    if(bool_vert_line){
      abline(v=f_center, col="blue")
    }
    dev.off()
  }
  
  # filename_data <- paste(prefix_filename, "_FFTdata_", spectrum_no_str, ".txt", sep="")
  write.table_with_header(format(cbind(frequency = freq_axis, psd = psd),
                                 digits = 10, nsmall=3, trim = TRUE),
                          file = filename_data,
                          header = "",
                          row.names = FALSE, quote = FALSE)
}


write.temperatures <- function(ps_save, filename_energies, dt, dt2 = 1e-4, number_of_Be_ions){
  
  # dt2 <- 1e-4
  # dt2 <- 1e-6
  R <- dt2/dt
  mydata <- ps_save[seq(from=2, to=nrow(ps_save), by=R),]
  mydata$T_p <- (qC2U0_p * mydata$zp^2 + 1/2 * m_p * mydata$vzp^2)/kb
  mydata$T_Be <- (qC2U0_Be * mydata$zBe^2 + 1/2 * m_Be * mydata$vzBe^2)/kb 
  
  mydata <- mydata[c("t","T_p", "T_Be")]
  
  # write.table_with_header(format(mydata , digits = 10, nsmall=3, trim = TRUE),
  #                         file = filename_energies,
  #                         header = "t\tT_p\tT_Be",
  #                         col.names = "t"
  #                         row.names = FALSE, quote = FALSE)
  
  write.table(format(mydata , digits = 10, nsmall=3, trim = TRUE),
                          file = filename_energies,
                          col.names = c("t", "T_p", "T_Be"),
                          row.names = FALSE, quote = FALSE)
  
}

# bla <- read.table("./2020-12-13_15h00_p_20Be_laser_detuning_scan_duty_cycle_1s/2020-12-13_15h00_p_20Be_laser_detuning_scan_duty_cycle_1s _energies 002 .txt",
#                   skip = 1)
# bla <- bla[which(bla[,1] > 1),]
# mean(bla[,3])/(20*kb)
# plot(bla[,1], bla[,3]/(20*kb), type="o", xlim=c(5,5.1))
# plot(bla[,1], bla[,2]/kb, type="l", xlim=c(5,30))

copy_sourcefiles <- function(complete_filename, filenameSubDir){
  
  if( (nchar(filenameSubDir)*2+90) > 255){
    stop("Filename too long!")
  }
  
  this_file_here <- complete_filename
  file.copy(this_file_here, paste("./",filenameSubDir, sep=""), overwrite=TRUE )
  match_pos <- regexpr("/.[^/]*\\.R", this_file_here)
  match_length <- attr(regexpr("/.[^/]*\\.R", this_file_here), "match.length")
  oldname <-substr(this_file_here, match_pos, match_pos + match_length)
  newname <- paste(substr(this_file_here, match_pos, match_pos + match_length -3), "_savecopy.R", sep="") 
  file.rename(from=paste("./", filenameSubDir, oldname, sep=""),
              to=paste("./", filenameSubDir, newname, sep=""))
  file.copy("./src/myRcppFile.cpp", paste("./",filenameSubDir, sep=""), overwrite=TRUE )
  file.rename(from=paste("./", filenameSubDir, "/myRcppFile.cpp", sep=""),
              to=paste("./", filenameSubDir, "/myRcppFile_savecopy.cpp", sep=""))
}

copy_resultfiles <- function(prefix_filename, filenameSubDir){
  file.copy("./monitorfile.txt", paste("./",filenameSubDir, sep=""), overwrite=TRUE )
  file.rename(from=paste("./", filenameSubDir, "/monitorfile.txt", sep=""),
              to=paste("./", filenameSubDir, "/monitorfile_savecopy.txt", sep=""))
  save.image(file=paste(prefix_filename, "_Rdata_savecopy.RData", sep=""))
}


# plot(Tdata_summary$Spectrum_number, Tdata_summary$T_Be/(10*Tdata_summary$T_p))




