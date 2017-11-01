


#========================================================================================
#  This file demonstrates the use of the program MultiSpect(), which
#  implements the method discussed in 
#  "Adaptive Bayesian Power Spectrum Analysis of Multivariate Nonstationary Time Series."
#
#   Content:
#   (1) Bivariate slowly-varying time series  (section 4.1)
#       (1a) Load the simulated data for bivariate slowly-varying time series
#       (1b) Run the Bayesian estimation process
#       (1c) Obtain and plot the number and location of partitions
#       (1d) Obtain and plot the surface of spectra and coherence
#       (1e) Obtain and plot the credible intervals for surfaces
#   (2) Multi-modal time series  (section 4.3)
#       (2a) Simulate data for bivariate slowly-varying time series
#       (2b) Run the Bayesian estimation process
#       (2c) Obtain and plot the number and location of partitions


#========================================================================================
### load packages and input functions
#========================================================================================
if("trust" %in% rownames(installed.packages()) == FALSE){install.packages("trust")} 
if("mvtnorm" %in% rownames(installed.packages()) == FALSE){install.packages("mvtnorm")}  
if("pracma" %in% rownames(installed.packages()) == FALSE){install.packages("pracma")} 
if("pscl" %in% rownames(installed.packages()) == FALSE){install.packages("pscl")}  
if("fields" %in% rownames(installed.packages()) == FALSE){install.packages("fields")} 
library(trust)
library(mvtnorm)
library(pracma)
library(pscl)
library(fields)
setwd("C:/Users/zeda/Dropbox/JASA submission/MultiSpect_R")
source("Beta_derive1.R");source("Beta_derive2.R")
source("birth.R"); source("death.R"); source("within.R")
source("chol_ind.R"); source("lin_basis_func.R")
source("Gradient1.R"); source("Gradient2.R")
source("Hamilt1.R"); source("Hamilt2.R")
source("postBeta1.R"); source("postBeta2.R")
source("whittle_like.R")
source("MultiSpec.R")
source("MultiSpect_surface.R")
source("MultiSpect_interval.R")
source("MultiSpect_partition.R")


### (1) 

#==============================================================================
#  (1a) Load the simulated data for bivariate slowly-varying time series
#==============================================================================

zt1 <- as.matrix(read.table("SlowVarying.txt"))
nobs <- dim(zt1)[1]
dim <- dim(zt1)[2]

#==============================================================================
#  (1b) Run the Bayesian estimation process
#==============================================================================
set.seed(20170531)
output1 <- MultiSpect(zt1, nloop=10000, nwarmup=2000, nexp_max=10, tmin=60, 
                     prob_mml=0.8, nbasis=10, tau_up_limit=10^4, sigmasqalpha=10^5,
                     init=3, nfreq=50, ee=0.1)

#==============================================================================
#  (1c) Obtain and plot the number and location of partitions
#==============================================================================
pos_prob <- MultiSpect_partition(zt1,output1)
pos_prob
#==============================================================================
#  (1d) Obtain and plot the surface of spectra and coherence
#==============================================================================
surface <- MultiSpect_surface(zt1,output1)
nfreq=50
par(mfrow=c(1,3))
image(1:nobs, (0:nfreq)/(2*nfreq), t(Re(surface$spectra1)),col=tim.colors(),zlim=c(0,8),xlab="Time", ylab="Frequency")
image(1:nobs, (0:nfreq)/(2*nfreq), t(Re(surface$spectra2)),col=tim.colors(),zlim=c(0,8),xlab="Time", ylab="Frequency")
image(1:nobs, (0:nfreq)/(2*nfreq), t(Re(surface$coh)),col=tim.colors(),zlim=c(0,1),xlab="Time", ylab="Frequency")

#============================================================================
# (1e) Obtain and plot the surface of spectra and coherence with credible intervals
#============================================================================
alphalevel <- 0.05
surface_interval <- MultiSpect_interval(zt1,output1,alphalevel)

par(mfrow=c(1,3))
image(1:nobs, (0:nfreq)/(2*nfreq), t(Re(surface_interval$spect11[,,2])),col=tim.colors(),zlim=c(0,10),main="lower bound",xlab="Time", ylab="Frequency")
image(1:nobs, (0:nfreq)/(2*nfreq), t(Re(surface_interval$spect11[,,1])),col=tim.colors(),zlim=c(0,10),main="mean",xlab="Time", ylab="Frequency")
image(1:nobs, (0:nfreq)/(2*nfreq), t(Re(surface_interval$spect11[,,3])),col=tim.colors(),zlim=c(0,10),main="upper bound",xlab="Time", ylab="Frequency")

image(1:nobs, (0:nfreq)/(2*nfreq), t(Re(surface_interval$spect22[,,2])),col=tim.colors(),zlim=c(0,12),main="lower bound",xlab="Time", ylab="Frequency")
image(1:nobs, (0:nfreq)/(2*nfreq), t(Re(surface_interval$spect22[,,1])),col=tim.colors(),zlim=c(0,12),main="mean",xlab="Time", ylab="Frequency")
image(1:nobs, (0:nfreq)/(2*nfreq), t(Re(surface_interval$spect22[,,3])),col=tim.colors(),zlim=c(0,12),main="upper bound",xlab="Time", ylab="Frequency")

image(1:nobs, (0:nfreq)/(2*nfreq), t(Re(surface_interval$coh21[,,2])),col=tim.colors(),zlim=c(0,1),main="lower bound",xlab="Time", ylab="Frequency")
image(1:nobs, (0:nfreq)/(2*nfreq), t(Re(surface_interval$coh21[,,1])),col=tim.colors(),zlim=c(0,1),main="mean",xlab="Time", ylab="Frequency")
image(1:nobs, (0:nfreq)/(2*nfreq), t(Re(surface_interval$coh21[,,3])),col=tim.colors(),zlim=c(0,1),main="upper bound",xlab="Time", ylab="Frequency")
par(mfrow=c(1,1))


### (2) 
#==============================================================================
#  (2a) Load the simulated data for bivariate slowly-varying time series
#==============================================================================
zt2 <- as.matrix(read.table("multimodal.txt"))
nobs <- dim(zt2)[1]
dim <- dim(zt2)[2]


#==============================================================================
#  (2b) Run the Bayesian estimation process
#==============================================================================
set.seed(20170531)
output2<- MultiSpect(zt2, nloop=10000, nwarmup=2000, nexp_max=10, tmin=100, 
                     prob_mml=0.8, nbasis=10, tau_up_limit=10^4, sigmasqalpha=10^5,
                     init=3, nfreq=50, ee=0.1)

#==============================================================================
#  (2c) Obtain and plot the number and location of partitions
#==============================================================================
pos_prob <- MultiSpect_partition(zt2,output2)
pos_prob

