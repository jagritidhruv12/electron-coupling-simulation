# Loading all the libraries we need
library(minpack.lm)
library("Hmisc")
library(Rcpp)
library(MySimulationsF)
library(foreach)
library(doParallel)
library(beepr)
library(zoo)

# Set working directory to the directory this file is located in
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# This loads another R-file where some useful functions are programmed
# -> serves to make this main program shorter
source("./utility_coupling.R")
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")


gc() # garbage collector, good but not necessary to do this at start

# Start time, gets printed in the RStudio console (or whatever you use)
print(paste("Start:", format(Sys.time(), "%d %b  %X")))
start_time <- proc.time()[[3]]


# Some constants for convenient calculations.
# !!! IMPORTANT !!!
# These are NOT used for the numerical calculations, this is only convenience for quick checks.
# The ones used in the numerical calculation are in the Cpp file or get changed by the user 
# further below, e.g. myrun$Dz_e2 <- ...
m_e1 <- 9.1093837139e-31
m_e2 <- 9.1093837139e-31
Lmult <- 1
myL <- Lmult*3e-3
f_det <- 626951.0
q <- 1.602176634e-19
Q <- 20e3
Dz_e1 <- 3.2e-3
Dz_e2 <- 3.2e-3
T_det <- 4.0
Tz_e1 <- 4
R_theo <- 2*pi*f_det*myL*Q
R <- 236e15
L <- 3e-3
C <- 2.14808329e-11
C_common <- 20e-12
qC2U0_e1 <- 1.29776307354e-14
qC2U0_e2 <-  1.29776307354e-14
kb <- 1.380649e-23
z0 <- sqrt(1.0 / 2 * kb * Tz_e1 / (qC2U0_e1))
w_z <- sqrt(2*qC2U0_e1/mp)
c <- 299792458.0
f0_e2Excit <- 957.39680e12
#photon_momentum <- 2.116056688185771e-27
# tau_e2_excitation <- 8.1e-9
tau_e2_excitation <- 8.1e-20
B0 <- 1.8
C2_e1 <- 18510
C2_e2 <- 18510
D4_e1 <- 1.283e9
D4_e2 <- 1.283e9
wc_e1 <- q/m_e1 * B0
wz_e1 <- sqrt(2*qC2U0_e1/m_e1)
wplus_e1 <-  wc_e1/2 + sqrt( wc_e1^2/4 - wz_e1^2/2 )
wminus_e1 <-   wc_e1/2 - sqrt( wc_e1^2/4 - wz_e1^2/2 )
fplus_e1 <- 1/(2*pi) * wplus_e1
fminus_e1 <- 1/(2*pi) * wminus_e1
wc_e2 <- q/m_e2 * B0
wz_e2 <- sqrt(2*qC2U0_e2/m_e2)
wplus_e2 <-  wc_e2/2 + sqrt( wc_e2^2/4 - wz_e2^2/2 )
wminus_e2 <-   wc_e2/2 - sqrt( wc_e2^2/4 - wz_e2^2/2 )
fplus_e2 <- 1/(2*pi) * wplus_e2
fminus_e2 <- 1/(2*pi) * wminus_e2
fz_e2 <- 1/(2*pi) * sqrt(2*qC2U0_e2/m_e2)



# change parameter to scan below here ------------------------------------------------------------


# This is the subdirectory where the data and a save-copy of this file will be stored 
filenameSubDir <- "2026-04-18_11h42_tutorial_simulation1" 


prefix_filename <- paste("./", filenameSubDir, "/", filenameSubDir, sep="")

# You are not allowed to put two simulations in the same folder.
if(dir.exists(file.path(this.dir, filenameSubDir)) ){
  stop("File directory already exists, please change directory name! The simulation has not been started.")
}else{
  dir.create(file.path(this.dir, filenameSubDir) )
}

this_file_here <- parent.frame(2)$ofile
copy_sourcefiles(this_file_here, filenameSubDir) # Copies a save file


# Important parameters. dt is the time step "\Delta t" and multiplier_billion_steps is how many steps you 
# want to simulate, here e.g. 1.2 billion, so in total 4e-9*6.0e9 s = 24s of real time.
# Note that number_of_steps/save_every_xth must be an integer number, so multiplier_billion_steps must be 
# a multiple of save_every_xth.
# e2 aware the sampling rates must be high enough so that the frequency you want to simulate is covered.
dt <- 4e-9
multiplier_billion_steps <- 0.06
number_of_steps <- 1e9*multiplier_billion_steps
save_every_xth <- 375


# This is the maximum number of cores you allow. My computer's process has 16 virtual cores, so I allow 15 
# parallel R threads. For this tutorial and to not crash your computer, I will limit it to 5
number_of_cores <- 3

# Estimation of how many Gigabyte this simulation needs (8 byte per double), 7 columns
print(paste("Estimated RAM needed:",
            format(8*number_of_steps/save_every_xth * 1e-9 * number_of_cores * 7, digits=2) , 
            "GB (not accurate for small numbers below 1GB)"))


# These are the parameters we want to go through.
p_absdetuning_arr <- c(-100, -50, -20) # First array: a proton tuned by 100, 50, and 20Hz to the left of the resonator
dnuz_e2 <- c(10, 50, 100) # Second array: Dip width of our beryllium ion

# Now we form our two arrays into a grid-matrix. The entries of this grid-matrix will later be used to specify
# the parameters for each run.
mygrid <- expand.grid(p_absdetuning_arr, dnuz_e2)


# this just does some checks, it is nothing that should (regularly) be changed
paramToScan <- mygrid[,1]
doSomeChecks(dt, number_of_steps, save_every_xth, FALSE, fplus_e2)

# stop("test")

# this sets up the parallel computing, it should not be  (regularly) changed
cl <- setupParallelComputing(number_of_steps, save_every_xth, paramToScan, number_of_cores)
gc()

# Try catch to also stop the cluster if an error occurred (curly brackets for more than one expression)
Tdata_summary <- tryCatch({
  Tdata_summary <- foreach(i=1:length(paramToScan), .combine=rbind) %dopar% {
    
    # i <- 1
    
    myrun <- new(SimulationRun) # initialize the simulation run object
    spectrum_no_str <- get_spectrum_no_str(i) # so that the first FFT is not called FFT1 but FFT001 (sorting)
    
    myrun$number_of_steps <- number_of_steps # number of loop iterations
    myrun$dt <- dt # time step
    
    myrun$save_every_xth <- save_every_xth # only every e.g. 375th step is saved
    myrun$bool_simulate_e1 <- TRUE # whether p shall be simulated at all
    myrun$bool_simulate_e2 <- TRUE # whether e2 shall be simulated at all
    myrun$bool_calculate_single_e2 <- TRUE # TRUE = independent e2 ions (computationally slow), FALSE = one effective ion
    myrun$number_of_e2_ions <- 1 # the number of independent e2 ions you want to simulate
 
    myrun$threadNo <- i # threadNo = loop iteration number
    
    
    myrun$T_init_e1 <- 4 # the initial electron 1 temperature
    myrun$T_init_e2 <- 4 # the initial electron 2 temperature 
 
    myrun$p_absdetuning <- 0 + mygrid[i,1] # the relative detuning of the p to the resonator,
    # this was the first column of our grid
    myrun$e2_absdetuning <- 0 # the relative detuning of the e2 to the resonator
    
    
    myrun$L <- 3.93e-3 # inductance of the coil
    myrun$f_det <- 345242 # resonator frequency
    Q <- 345242/30 # resonator Q-factor, which is used to calculate R
    myrun$R <- 2*pi*myrun$f_det * myrun$L * Q # calculation of R
  
    dnuz_p <- 9 # Proton shall have 9Hz here
    myrun$Dz_e1 <- sqrt(1/(2*pi) * myrun$R/m_e1 * q^2/dnuz_p)
  
    myrun$Dz_e2 <- sqrt(1/(2*pi) * myrun$R/m_e2 * q^2/mygrid[i,2] * myrun$number_of_e2_ions) 
  
    
    
    # myrun$stability_voltage_p # voltage stability
    myrun$Dz_cap_p # capacitive trap size
  
    # calculate some stuff based on parameters
    qC2U0_e1 <- 2 * pi^2 * m_e1 * (myrun$f_det + myrun$p_absdetuning)^2
    qC2U0_e2 <- 2 * pi^2 * m_e2 * (myrun$f_det + myrun$e2_absdetuning)^2
    U0_e1 <- qC2U0_e1/(q*C2_e1)
    U0_e2 <- qC2U0_e2/(q*C2_e2)
    
    fz_p <- myrun$f_det + myrun$p_absdetuning
    fz_e2 <- myrun$f_det + myrun$e2_absdetuning
    
    number_of_e2_ions <- myrun$number_of_e2_ions
    bool_calculate_single_e2 <- myrun$bool_calculate_single_e2
    
    C_common <- myrun$C_common
    
    
    
    ps_save <- myrun$performLoop() # This calls the performLoop()-function in the Cpp-program, which,
    # as the name suggests, performs the computationally expensive loop. ps for R space 
    
    # Append the time column (doing it here saves a little RAM)
    ps_save <- cbind(seq(from=0, to=(number_of_steps-1)*dt, by=save_every_xth*dt), ps_save)
    
    C <- myrun$C # just to be sure we have the correct values, we can call them again from the myrun-object
    L <- myrun$L
    
    remove(myrun) # this removes the myrun-object and "finishes" the Cpp-part
    gc()
    
    ps_save <- as.data.frame(ps_save) # just some data conversion
    
    # IMPORTANT! If you change the variables that are fed back from the Cpp-program, you must change 
    # the ps_save-object accordingly, including the line below (otherwise everything will be mixed up). 
    # Time, Electron 1 z, Electron 1 v_z, Resonator voltage, Current through coil, e2 cloud total energy
    colnames(ps_save) <- c("t", "zp", "vzp", "Ures", "IL", "ze2", "vze2", "Ee2")  
    
    
    ps_save$T_e2 <- ps_save$Ee2/kb
    
    
    # The rest of the loop just calculates the FFT spectra and saves the data & plots. I will not comment
    # every detail here, because it should be straightforward
    dt2 <- 1e-2
    K <- floor(dt2/(dt*save_every_xth))
    Tdata <- ps_save
    Tdata$T_e1 <- (qC2U0_e1 * Tdata$zp^2 + 1/2 * m_e1 * Tdata$vzp^2)/kb
    Tdata$T_e2 <- Tdata$Ee2 / kb 
    if(number_of_e2_ions > 0){
      if(!bool_calculate_single_e2){
        Tdata$T_e2 <- Tdata$T_e2*number_of_e2_ions
      }
    }else{
      Tdata$T_e2 <- 0
    }
    Tdata$T_RLC <- (1/2 * C * Tdata$Ures^2 + 1/2 * L * Tdata$IL^2)/kb
    Tdata$T_total <- Tdata$T_e1 + Tdata$T_e2 
    
    
    Tcolnames <- c("t","T_e1", "T_e2", "T_RLC",  "T_total")
    Tdata <- Tdata[Tcolnames]
    
    
    plottext <- paste("i = ", i, "\n mygrid[i,1] = ", mygrid[i,1], ", mygrid[i,2] = ", mygrid[i,2])
    
    
    png(paste(prefix_filename, "_energy_exchange", spectrum_no_str, ".png", sep=""),
        width=1080*2, height=720*2, pointsize=20*2)
    plot(Tdata$t, Tdata$T_total, type="l", ylim=c(0, max(Tdata$T_total)*1.1), xlab = "Time [s]", ylab="Energies / kB (K)",
         col="green")
    lines(Tdata$t, Tdata$T_e2, type="l", col="blue")
    lines(Tdata$t, Tdata$T_e1, type="l")
    lines(Tdata$t, Tdata$T_RLC, type="l", col="red")
    lines(Tdata$t, Tdata$T_com, type="l", col="purple")
    legend("topleft", c("Electron 1", "Electron 2", "Resonator", "C_com", "Sum"),
           col=c("black", "blue", "red", "purple", "green"), lty=c(1,1,1,1) )
    mtext(plottext)
    dev.off()
    
    png(paste(prefix_filename, "_energy_exchangeA", spectrum_no_str, ".png", sep=""),
        width=1080*2, height=720*2, pointsize=20*2)
    plot(Tdata$t, Tdata$T_e2/number_of_e2_ions, type="l", ylim=c(0, max(c(Tdata$T_e1, Tdata$T_e2/number_of_e2_ions))*1.1),
         xlab = "Time [s]", ylab="Energies / kB (K)", col="blue")
    lines(Tdata$t, Tdata$T_e1, type="l", col="black")
    legend("topright", c("Electron 1", "Electron 2"), col=c("black", "blue"), lty=c(1,1) )
    mtext(plottext)
    dev.off()
    
    png(paste(prefix_filename, "_energy_exchangeB", spectrum_no_str, ".png", sep=""),
        width=1080*2, height=720*2, pointsize=20*2)
    plot(Tdata$t, Tdata$T_e2, type="l", ylim=c(0, max(c(Tdata$T_e1, Tdata$T_e2))*1.02),
         xlab = "Time [s]", ylab="Energies / kB (K)", col="blue")
    lines(Tdata$t, Tdata$T_e1, type="l", col="black")
    legend("topright", c("Electron 1", "Electron 2"), col=c("black", "blue"), lty=c(1,1) )
    mtext(plottext)
    
    dev.off()
    
    
    
    filename_energies <- paste(prefix_filename, "_energies", spectrum_no_str, ".txt", sep="")
    
    write.table(format(Tdata , digits = 6, nsmall=3, trim = TRUE),
                file = filename_energies,
                col.names = Tcolnames,
                row.names = FALSE, quote = FALSE)
    
    
    
    
    plottext <- paste("i = ", i, "\n mygrid[i,1] = ", mygrid[i,1], ", mygrid[i,2] = ", mygrid[i,2])
    myheader <- ""
    
    
    
    fft_range_plus <- 400
    fft_range_minus <- 400
    ##################################################
    # ----------------
    #################
    
    plotfilename <- paste(prefix_filename,"_FFT_", spectrum_no_str ,".png", sep="")
    plotfilename2 <- paste(prefix_filename,"_FFTnarrow_", spectrum_no_str ,".png", sep="")
    filename_data <- paste(prefix_filename, "_FFTdata_", spectrum_no_str, ".txt", sep="")
    
    plot_and_save_FFT(ps_save$Ures, dt, save_every_xth, fz_p, 200, 10,
                      plotfilename, spectrum_no_str, plottext, filename_data)
    
    
    ## ----
    
    plotfilename <- paste(prefix_filename,"_FFTzp_", spectrum_no_str ,".png", sep="")
    filename_data <- paste(prefix_filename, "_FFTzp_data_", spectrum_no_str, ".txt", sep="")
    
    plot_and_save_FFT(ps_save$zp, dt, save_every_xth, fz_p, 400, 10,
                      plotfilename, spectrum_no_str, plottext, filename_data)
    
    # -----
    plotfilename <- paste(prefix_filename,"_FFTze2_", spectrum_no_str ,".png", sep="")
    filename_data <- paste(prefix_filename, "_FFTze2_data_", spectrum_no_str, ".txt", sep="")
    
    plot_and_save_FFT(ps_save$ze2, dt, save_every_xth, fz_e2, 400, 10,
                      plotfilename, spectrum_no_str, plottext, filename_data)
    
    
    
    
    myrow <- cbind(i, mean(Tdata$T_e1), mean(Tdata$T_e2), mean(Tdata$T_RLC) )
    
    rm(ps_save, Tdata)
    gc()
    
    
    myrow
    
    
    
  } # bracket for foreach()-loop
  Tdata_summary
}, # bracket for tryCatch {}-wrapper
finally={
  #stop cluster. Try catch to also stop cluster if an error occurred
  stopCluster(cl)
  
  
  gc()
}
) # bracket for tryCatch()-function    




# The rest of the script merely saves and plots some summary data 
if(ncol(Tdata_summary) > 4){
  colnames(Tdata_summary) <- c("Spectrum_number", "T_e1", "T_e2", "T_res", "sd_U0_e1", "sd_U0_e2", "rel_diff_U0")
}else{
  colnames(Tdata_summary) <- c("Spectrum_number", "T_e1", "T_e2", "T_res")
}
Tdata_summary <- as.data.frame(Tdata_summary)




# End time
print(paste("End:", format(Sys.time(), "%a %d %X")))
end_time <- proc.time()[[3]]
print(paste("Seconds passed:", end_time-start_time))
print("Finished loop!")

filename2 <- paste(prefix_filename, "_Tdata_summary.txt",sep="")
write.table(Tdata_summary, filename2, row.names = FALSE, quote = FALSE)



plot(Tdata_summary$Spectrum_number, Tdata_summary$T_res, xlab="Spectrum number", ylab="T resonator (K)",  pch=16,
     ylim=range(c(Tdata_summary$T_res, Tdata_summary$T_e1, Tdata_summary$T_e2)), type="o")
lines(Tdata_summary$Spectrum_number, Tdata_summary$T_e1, type="o", pch=16, col="red")
lines(Tdata_summary$Spectrum_number, Tdata_summary$T_e2, type="o",  pch=16, col="blue")
legend("topright", c("Resonator", "Electron 1", "Electron 2"), col=c("Black", "Red", "Blue"), lty=c(1,1,1))



png(paste(prefix_filename,"_summary_temperatures.png", sep=""), width=1080*2, height=720*2, pointsize=20*2)
plot(Tdata_summary$Spectrum_number, Tdata_summary$T_res, xlab="Spectrum number", ylab="T resonator (K)",  pch=16,
     ylim=range(c(Tdata_summary[,-1])), type="o")
lines(Tdata_summary$Spectrum_number, Tdata_summary$T_e1, type="o", pch=16, col="red")
lines(Tdata_summary$Spectrum_number, Tdata_summary$T_e2, type="o",  pch=16, col="blue")
legend("topright", c("Resonator", "Electron 1", "Electron 2"), col=c("Black", "Red", "Blue"), lty=c(1,1,1))
dev.off()



png(paste(prefix_filename,"_summary_T_e1.png", sep=""), width=1080*2, height=720*2, pointsize=20*2)
plot(Tdata_summary$Spectrum_number, Tdata_summary$T_e1, pch=16, xlab="Spectrum number", ylab="T proton (K)")
dev.off()

png(paste(prefix_filename,"_summary_T_e2.png", sep=""),
    width=1080*2, height=720*2, pointsize=20*2)
plot(Tdata_summary$Spectrum_number, Tdata_summary$T_e2, pch=16, xlab="Spectrum number", ylab="T e2 (K)")
dev.off()

png(paste(prefix_filename,"_summary_T_res.png", sep=""),
    width=1080*2, height=720*2, pointsize=20*2)
plot(Tdata_summary$Spectrum_number, Tdata_summary$T_res, pch=16, xlab="Spectrum number", ylab="T resonator (K)")
dev.off()


copy_resultfiles(prefix_filename, filenameSubDir) # copy the monitorfile in the folder

beep(sound=1) # A microwave beep sound so you know your simulations are done ;) 


