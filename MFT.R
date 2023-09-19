
mft <- function(series, tol=1.e-3) {
  MeanSeries <- mean(series)
  series <- series - MeanSeries
  times <- seq(start(series)[1], end(series)[1],  1/frequency(series))
  dt = 1/frequency(series)
  N <- length(series)
  freqs <- 2*pi/(N*dt)*seq(0, floor(N/2) - 1)
  deltaf <- (2*pi)/(N*dt)

  ic = complex(1,0,1)

  found_frequency <- c()
  found_amplitude <- c()
  found_phases <- c()
  var0 <- var(series)
  tempseries <- series
  
  
  while (var(tempseries)/var(series) > tol) {
    residual_variance <- var(tempseries)/var(series)
    print ('iteration')
    print (residual_variance)
    Fseries <- Mod(fft(tempseries)[1:floor(N/2)]/N)
    isigper <- which.max(Fseries)
    Fmax   <- freqs[  isigper ]
    finetune <- seq(Fmax-deltaf/2., Fmax+deltaf/2., length=201)
    # ici il faut reflechire a la normalisation. Sinon il selectionne une mauvaise frequenc. 
     sumcos2 <- sapply(finetune, function(f) sum(cos(f*times)^2))
     sumsin2 <- sapply(finetune, function(f) sum(sin(f*times)^2))
     sumsincos <- sapply(finetune, function(f) sum(sin(f*times)*cos(f*times)))
     alpha <- sapply(finetune, function(f) sum(cos(f*times)%*%tempseries))
     beta <- sapply(finetune, function(f) sum(sin(f*times)%*%tempseries))
  
     A <- rep(0, length(finetune))
     B <- rep(0, length(finetune))
     squareamps <- rep(0, length(finetune))
     phases  <- rep(0, length(finetune))
  
      for (j in seq(along=finetune)){
        MAT <- matrix(c(sumcos2[j], sumsincos[j], sumsincos[j], sumsin2[j]), 2,2)
#         if ((j > 99) && (j < 103))  {print(MAT)}
       AB <- solve(MAT,c(alpha[j],beta[j])) 
       A[j] <- AB[1]
       B[j] <- AB[2]
       squareamps[j] <- AB[1]^2 + AB[2]^2
     }
  
#     plot(finetune, squareamps, type='l', xlim=c(0.78, 0.80), ylim=c(0.30, 0.35))
#     abline (v=iomega)
#     abline (h=a^2+b^2)
  
    argmax <- which.max(squareamps)
    fmax <- finetune[argmax]
    amp  <- sqrt(squareamps[argmax])
    phase <- Arg(A[argmax] - ic*B[argmax])
  
    sig_to_subtract <- ts(
         A[argmax] * cos(fmax*times) +
         B[argmax] * sin(fmax*times), deltat=dt)
  
    found_frequency <- c(found_frequency, fmax)
    found_amplitude <- c(found_amplitude, amp)
    found_phases    <- c(found_phases, phase)
    tempseries <- tempseries - sig_to_subtract
#     plot(tempseries)
  }
  
  found_periods <- 2*pi / found_frequency
  
  
  O2 <- order(found_amplitude, decreasing=TRUE)
  
  found_phases <- found_phases[O2]
  found_amplitude <- found_amplitude[O2]
  found_periods <- found_periods[O2]
  found_frequency <- found_frequency[O2]
  
  return(list(A=found_amplitude, O=found_frequency, P=found_phases, Period=found_periods))
}
  



