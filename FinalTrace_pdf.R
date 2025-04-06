# Calculate the pdf of distance traveled
# Assume a power law velocity pdf
#  distribution parameters
#    vmax - maximum velocity (m/s)
#    N    - dimensionless shape parameter

library(DescTools)

# powerlaw velocity distribution functions
pl.vpdf <- function(v, vm, N) { # power law pdf
  return((N/vm)*(v/vm)^(N-1)) # (s/m)
}
pl.vcdf<- function(v, vm, N) { # power law cdf
  return((v/vm)^N) # (dimensionless)
}
pl.vinvcdf<- function(P, vm, N) { # power law inverse cdf
  return(vm*P^(1/N)) # (m/s)
}

pl.Nest <- function(vavg, vm) { # estimate N from average and max v
  return(vavg/(vm-vavg)) # (dimensionless)
  }

# Power law velocity statistics
pl.vmean <- function(vm, N) { # power law mean velocity
  return(vm*N/(N+1)) # (m/s)
}
pl.vvar  <- function(vm, N) { # power law variance
  return(vm^2 *N/((N+1)^2 * (N+2))) # (m2/s2)
}

# power law distance statistics
pl.xpdf <- function(x, xm, N) { # power law pdf
  return((N/xm)*(x/xm)^(N-1)) # (m)
}
pl.xcdf<- function(v, vm, N) { # power law pdf
  return((v/vm)^N) # (dimensionless)
}

# power law time/slowness statistics
pl.smean <- function(vm,N) { # power law mean slowness s=1/v
  smin <- 1/vm # (w/m)
  # initialize so dim of savg is the sme as dim of N
  savg <- rep(Inf, length(N))
  Ng1 <- (N>1)
  savg[Ng1] <- smin*N[Ng1]/(N[Ng1]-1)
  return(savg) # (s/m)
}

# power law random generator
pl.rv <- function(k, vm, N) {
  return(pl.vinvcdf(runif(k), vm, N))
}

# cloud functions
clouds.cdf.calc <- function(trand, ds=1) {
  # cloud cdf calculation downsampled to ns points
  #    trand has times of travel in a new col for each added mixing
  num <- length(trand[,1]) # number of observations
  tot <- apply(trand, 2, sort) # sort columns of trand
  probs <- ((1:num)-0.5)/num # corresponding probabilities
  # downsize 
  nkeep <- seq(1,num,by=ds) # subscripts to keep
  tot <- tot[nkeep,]
  probs <- probs[nkeep]
  clouds.cdf <- data.frame(tot, probs)
  #names(clouds.cdf <- c())
  return(clouds.cdf)
} # end clouds.cdf.calc

clouds.pdf.calc <- function(cdf, wingwidth=2) {
  # calculate the pdf from the derivative of the cdf
  # cdf is a data frame with time-of-travel (tot) in the first m cols
  #   and probability in final column
  # wingwidth determines the width of the window used to 
  #   calculate the derivative of P with respect to tot, 
  
  n <- dim(cdf)[1]   # number of observations in each cloud
  m <- dim(cdf)[2]-1 # number of tot clouds
  P <- cdf$probs     # probabilities
  w <- (2*wingwidth)+1 # width of sliding window
  dP <- rep(NA, n)   # initialize dp
  dP[(wingwidth+1):(n-wingwidth)] <-  # delta P
    P[w:n] - P[1:(n-w+1)]
  pdf <- cdf # initialize pdf 
  pdf[1:wingwidth, 1:m] <- NA       # derivatives near ends are NA
  pdf[(n-wingwidth+1):n, 1:m] <- NA
  for (j in 1:m) {
    for (i in (wingwidth+1):(n-wingwidth)) {
      dtot <- cdf[i+wingwidth,j]-cdf[i-wingwidth,j]
      pdf[i,j] <- # delta P / delta tot
        dP[i] / dtot
    } # end for i
  } # end for j
  return(pdf)
} # end clouds.pdf.calc

cloud.pdf.smooth <- function(pdf, cdf, wingwidth=2, knots=20) {
  # smmooth the pdf functions in pdf dataframe
  n <- dim(pdf)[1]   # number of observations in each pdf
  m <- dim(pdf)[2]-1 # number of pdf distributions
  w <- wingwidth     # needed to remove NA from smoothing
  r <- (w+1):(n-w)   # index range for smoothing
  pdf.sm <- pdf
  for (j in 1:m) { # loop through the distributions
    smodel <- 
      smooth.spline(cdf[r,j], pdf[r,j], nknots=knots)
    pdf.sm[r,j] <- smodel$y
  } # end for j
  return(pdf.sm) # return the smoothed distributions
}

# UNDER CONSTRUCTION
# cloud.pdf.normalize <- function(pdf, )

# ----------------- Main ----------------------
# CONVENTION - dimensionless UPERCASE, dimensional lowercase

# STATE VARIABLES
#   t  - time of travel (T)
#   x  - distance downstream
#
# PARAMETERS
#   N  - power law shape factor (dimensionless)
#   vm - maximum velocity (M/T)
#   l  - mixing length (L)
#   xf - final x, sampling site (L)
#
# CALCULATED
#   vavg - power law average velocity (L/T)
#   tavg - power law avg time of travel over mixing length l (T)
#   tmin - power law min time of travel over mixing length l (T)
#
# NONDIMENSIONAL
#   V - velocity = v/vm
#   X - downstream distance = x/L
#   T - time of travel = t/tavg

# number of random values generated in Monte Carlo simulation
  num <- 100000

#     use l=1, vm=1 for dimensionless results
l  <- 1
vm <- 1

# user input
xf <- readline('Enter distance to sampling site (xf>0) = ')
N  <- readline('                    Enter shape factor = ')
xf <- as.numeric(xf)
N <- as.numeric(N)


# number of times to move including final partial move 
m <- ceiling(xf/l) 
# results of simulations are saved in trand
trand <- matrix(nrow = num, ncol = m)

for (j in 1:num) { # collect num results
  # initialize
  x <- 0
  t <- 0
  for (i in 1:m) {
    if ((x+l) < xf) { # this step will not pass the sampling site
      xstep <- l # take a full step
    } else { # final step
      xstep <- xf - x
    } # end if
    x <- x+xstep
    v <- pl.rv(1, vm, N) # new velocity for this step
    tstep <- xstep/v # time to take this step
    t <- t+tstep # total time to travel to this point
    trand[j,i] <- t # save jth time of travel at ith mixing
    } # end for i
  } # end for j
summary(trand)
tavgcalc <- xf*pl.smean(vm,N)
print(paste('calculated average time =', tavgcalc))
print(paste(' simulated average time =', mean(trand[,m])))
# initial plot
hist(trand[,m][trand[,m]<quantile(trand[,m],probs = 0.95)],
     freq = FALSE, xlab = 'Travel time', ylab = 'pdf',
     main = paste('Travel', xf, 'mixing lengths, N=',N))
# cloud calculation
  clouds.cdf <- apply(trand, 2, sort) # sort columns of trand
  probs <- ((1:num)-0.5)/num # corresponding probabilities
# locate the tails
  # initialize
  cloud.p05 <- rep(NA,m)
  cloud.p95 <- rep(NA,m)
  width      <- rep(NA,m)
for (i in 1:m) { # for each mixing
  cloud.p05[i] <- quantile(trand[,i], probs = 0.05)
  cloud.p95[i] <- quantile(trand[,i], probs = 0.95)
  width[i]      <- cloud.p95[i] - cloud.p05[i]
  # determine the peak
  
} # end for i
  plot(1:m,width, col='red', 
       xlab = 'Distance', ylab = 'Cloud width', 
       main = paste('Travel', xf, 'mixing lengths, N=',N))
  lines(1:m,width[m]*seq(1/m,1, by = 1/m), col='gray')
  lines(1:m,width[m]*(seq(1/m,1, by = 1/m)^(1/2)), col='green')
  # create cdf
  cdf <- clouds.cdf.calc(trand,ds=100)
  
# calculate the pdf
  pdf <- clouds.pdf.calc(cdf, wingwidth=2)
  pdf.sm <- cloud.pdf.smooth(pdf, cdf, wingwidth=2, knots=20)

  # plot the pdfs
  plot(cdf[,1], pdf.sm[,1],
       ylim = c(0, max(pdf.sm[,2]*1.1, na.rm = TRUE)),
       xlim = c(l, 2*xf), type = 'l')
  for (j in 2:m) {
    lines(cdf[,j], pdf.sm[,j])
  }
  
  