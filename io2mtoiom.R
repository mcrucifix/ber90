I2OM <- read.table('iom.dat')
colnames(I2OM) <- c('index','sigmaprime','Nprime','deltaprime')


# I2OM <- I2OM[-1,]  # suppress constant term 

# need to think more about this. 

# transform into radians / year and radians

I2OM$sigmaprime <- I2OM$sigmaprime/(60*60*360)*2*pi
I2OM$Nprime  <- I2OM$Nprime / 1.e8
I2OM$deltaprime <- I2OM$deltaprime * 2*pi / 360. 


n <- length(I2OM$Nprim)

# I will privilege here transparent writing, even if less efficient. Fortran compiler would be better

# line 1

TERM1 <- list()
TERM1$sigma <- I2OM$sigmaprime
TERM1$N <- I2OM$Nprime
TERM1$delta <- I2OM$deltaprime

for (i in seq(n)) TERM1$N[i] <- TERM1$N[i] * ( 2 - I2OM$Nprime[i]^2 - 2 * sum(I2OM$Nprime[-i]^2) ) 

# this corresponds to lines p7505810 in his code 
# line 2 (TERM2)

TERM2 <- list()

TERM2$N <- - outer(I2OM$Nprime, I2OM$Nprime, function(Ni,Nj) Nj^2*Ni);
TERM2$sigma <- outer(I2OM$sigmaprime, I2OM$sigmaprime, function(Si,Sj) (2*Sj-Si))
TERM2$delta <- outer(I2OM$deltaprime, I2OM$deltaprime, function(Si,Sj) (2*Sj-Si))

for (h in seq(ncol(TERM2$N))) TERM2$N[h,h] = NA
for (h in seq(ncol(TERM2$sigma))) TERM2$sigma[h,h] = NA
for (h in seq(ncol(TERM2$delta))) TERM2$delta[h,h] = NA

TERM2$N <- as.numeric(TERM2$N)
TERM2$sigma <- as.numeric(TERM2$sigma)
TERM2$delta <- as.numeric(TERM2$delta)


TERM2$N <- TERM2$N[which(!is.na(TERM2$N))]
TERM2$sigma <- TERM2$sigma[which(!is.na(TERM2$sigma))]
TERM2$delta <- TERM2$delta[which(!is.na(TERM2$delta))]


# line 2 (TERM3)

TERM3 <- list()
TERM3$N <- -2 * outer ( outer(I2OM$Nprime, I2OM$Nprime[-n]), I2OM$Nprime)

TERM3$sigma<- outer( outer(I2OM$sigma, I2OM$sigma[-n], function(si, sj) sj-si), I2OM$sigma, function(sj,sk) sj+sk)

TERM3$delta<- outer( outer(I2OM$delta, I2OM$sigma[-n], function(si, sj) sj-si), I2OM$delta, function(sj,sk) sj+sk)

for (i in seq(n)) {
  for (j in seq(n-1)) {
    for (k in seq(n)) {
      if ((k >= j) || (k == i) || (j == i) ){
        TERM3$sigma[i, j, k ] = NA;
        TERM3$delta[i, j, k ] = NA;
        TERM3$N[i, j, k ] = NA;
      }
    }}}

TERM3$sigma <- as.numeric(TERM3$sigma)
TERM3$delta <- as.numeric(TERM3$delta)
TERM3$N <- as.numeric(TERM3$N)


TERM3$N <- TERM3$N[which(!is.na(TERM3$N))]
TERM3$sigma <- TERM3$sigma[which(!is.na(TERM3$sigma))]
TERM3$delta <- TERM3$delta[which(!is.na(TERM3$delta))]


IOM <- list(
            N = c(TERM1$N, TERM2$N, TERM3$N), 
            sigma = c(TERM1$sigma, TERM2$sigma, TERM3$sigma), 
            delta = c(TERM1$delta, TERM2$delta, TERM3$delta))


# THE IOM LIST IS VALIDATED. AT LEAST IN AMPLITUDES, AND FOR THE FIRST SERIES OF 80 TERMS
tol <- 1.e-7

nkeep <- which (IOM$N > tol)

IOM <- list( N = IOM$N[nkeep], 
             sigma = IOM$sigma[nkeep], 
             delta = IOM$delta[nkeep]) 



#### development of obliquity

epsilonbar = 23.44579 * pi / 180
psibar     = 0.0002445400207  + 0.0002443015 - 0.0002445400 
zeta       = 1.964 * pi / 180. 

# thesis Berger p. 110 = Sharaf Budnikova

fi <-  psibar + IOM$sigma
ci <-  psibar / fi 
deltaprime <- IOM$delta + zeta 

tane = tan(epsilonbar)
cote = 1./tane

nf <- length(ci)

C1f <- -ci
C2fik <- -0.5 * psibar * outer (seq(nf), seq(nf), function(i,k) {
                       1./(fi[i]+fi[k])*((ci[i]^2 + ci[k]^2 - 2) * tane +
                              (ci[i]+ci[k]) * cote)  })


C2fii = -.25 *(ci * (ci^2 - 1) * tane + ci^2 * cote )

C2fikprime <- -0.5 * psibar * outer (seq(nf), seq(nf), function(i,k) {
                       1./(fi[i]-fi[k])*((ci[i]^2 - ci[k]^2 + 2*ci[k] - 2*ci[i]) * tane +
                              (ci[i]-ci[k]) * cote)  })


## ATTENTION CECI DOIT ETRE CORRIGE
## C'EST EN FONCTION DE CIF

## mais je dois revoir la these pasge 113
## et convertir les cif en ci. 

## pour ca il faut calculer les Di, qui ne sont pas encore calcules !!!! 


#### DEVELOPPEMENT FOR PSI         

# P0 et ell ont ete repris de sa these, sauf que pour ell il y a peut etre un facteur (1-e_0)3/2

P0 = 17.3919
ell = 54.9066 
EPI <- read.table('epi.dat')
colnames(EPI) <- c('index','g','M','beta')
EPI$g <- EPI$g/(60*60*360)*2*pi
EPI$M  <- EPI$M / 1.e8
EPI$beta <- EPI$beta * 2*pi / 360. 

ng <- length(EPI$g)




# precession




D1f <- ci * (cote + (ci - 1) * tane)

D2f <- 0.25 * ci * (ci^2 + ci - 1) + 0.125 * ci^2*(ci-1)^2*(tane^2) + .5 * ci^2 * cote^2 

# ligne p7507950

D2fik <- psibar * outer(fi,fi, function(fi, fk){1/(fi + fk)}) * (
          0.5 * outer (ci, ci,  function(ci,ck) {ci^2+ck^2+ci+ck-ci*ck-2} )
         + C2fik * tane
         - 0.5 *   outer (ci, ci,  function(ci,ck) {ci*(ci-1)+ck*(ck-1)}) * tane^2
         + outer (ci,ci, '+')*cote^2 )


D2fikprime <- psibar * outer(fi,fi, function(fi, fk){1/(fi - fk)}) * (
          0.5 * outer (ci, ci,  function(ci,ck) { 5 * (ci+ck) - (ci^2+ck^2) - ci*ck - 6} )
         + C2fikprime * tane
         + 0.5 *  outer (ci*(ci-1), ci*(ci-1),'+') * tane^2) 







# must define here P0, ell but that's another part of the code. 

Dfseconde = 3 * psibar * P0 / ell * outer(EPI$g, EPI$g,  function(gi,gk) { 1./(gi-gk)})


D1 <- D1f - ci
D2 <- D2f - (ci - 0.5)*cote^2 - 0.5 * (ci^2 - 0.5)
D2ik <- D2fik - (outer(ci, ci, '+')-1) * cote^2 - 0.5 *(outer(ci^2, ci^2,'+')-1)
D2ikprime <- D2fikprime - 0.5 * ( outer(ci^2, ci^2,'-')) + outer(ci, ci, '+') 

# psi still need to be coded. eq. 64 of his thesis p. 112


# computes the big Ci

# remember that we have the 
#### C1f
#### C2f 
#### C2fii 
#### C2fprime
 
C1 = C1f + 1
C2ii = C2fii + 0.5*D1f - 0.25 * cote
C2ik = C2fik + 0.5 * outer(D1f, D1f, '+') - 0.50 * cote
C2ikprime = C2fikprime - 0.5 * outer (D1f, D1f, '+') + 0.50 * cote


ETERM1 <- list(  N = C1 * IOM$N, sigma = fi, delta = deltaprime ) 

# ETERM1 est validé

ETERM2 <- list ( N = C2ii * IOM$N^2 , sigma = 2*fi, delta = 2 * deltaprime ) 

# ETREM2 est validé

ETERM3 <- list ( N = C2ik * outer(IOM$N, IOM$N), sigma = outer(fi, fi, '+'), delta = outer(deltaprime, deltaprime, '+'))

for (i in seq(length(IOM$N))) {
  for (k in seq(length(IOM$N))) {
     if ( i >= k ) {
       ETERM3$N[i,k] = NA; 
       ETERM3$sigma[i,k] = NA;
       ETERM3$delta[i,k] = NA}}}


ETERM3$N <- ETERM3$N[which(!is.na(ETERM3$N))]
ETERM3$sigma <- ETERM3$sigma[which(!is.na(ETERM3$sigma))]
ETERM3$delta <- ETERM3$delta[which(!is.na(ETERM3$delta))]

# regarder ETERM4 : il est manifestement faux
# il y a aussi un mysterieux facteur 0.95, à moins que ce soit
# un facteur 2 sur omega i / vs i2

ETERM4 <- list ( N = C2ikprime * outer(IOM$N, IOM$N), sigma = outer(fi, fi, '-'), delta = outer(deltaprime, deltaprime, '-'))

for (i in seq(length(IOM$N))) {
  for (k in seq(length(IOM$N))) {
     if ( i >= k ) {
       ETERM4$N[i,k] = NA; 
       ETERM4$sigma[i,k] = NA;
       ETERM4$delta[i,k] = NA}}}


ETERM4$N <- ETERM4$N[which(!is.na(ETERM4$N))]
ETERM4$sigma <- ETERM4$sigma[which(!is.na(ETERM4$sigma))]
ETERM4$delta <- ETERM4$delta[which(!is.na(ETERM4$delta))]




OBLIQUITY <- list ( N = c(ETERM1$N, ETERM2$N, ETERM3$N, ETERM4$N), 
                    sigma = c(ETERM1$sigma, ETERM2$sigma, ETERM3$sigma, ETERM4$sigma), 
                    delta = c(ETERM1$delta, ETERM2$delta, ETERM3$delta, ETERM4$delta))

# keep the first 200 terms

# we could also transform the negative $N$ into positive by changing the phase (all cosines)

# again one step at a time. ... .

nkeep = 2000
OR <- order(abs(OBLIQUITY$N), decreasing=TRUE)[seq(nkeep)]
OBLIQUITY <- with(OBLIQUITY, list(N = N[OR], sigma=sigma[OR], delta=delta[OR]))

# with this, there seems be no duplicates (but maybe there will be once we have recitfied the N's)
# however we have a problem becase the sigma0 generates periods = psibar

which(duplicated(OBLIQUITY$sigma))

2*pi/ OBLIQUITY$sigma[seq(100)]

require(palinsol)
data(BER90)
par(mfrow=c(1,2))

plot(abs(OBLIQUITY$sigma)/(2*pi), OBLIQUITY$N, type='n', ylim=c(-1e-2, 0), xlim=c(0e-5, 4e-5))
lines(ETERM1$sigma/(2*pi), ETERM1$N, col='red', type='h')
lines(ETERM2$sigma/(2*pi), ETERM2$N, col='blue' , type='h')
lines(ETERM3$sigma/(2*pi), ETERM3$N, col='gray', type='h')
lines(ETERM4$sigma/(2*pi), ETERM4$N, col='green', type='h')
grid()
with(BER90$Table1,  plot(Rate/360/3600 , Amp*2*pi/60/60/360, type='h', col='black', ylim=c(-1e-2, 0), xlim=c(0e-5, 4e-5)))
grid()

# ok this is more or less validated. 

#plot(ETERM1$sigma/(2*pi), C1f, type='h')

# attention at this point we haven't tried to reorganise and chase doublons. 
# there will be many of those, but one step at a time. 
# make all frequencies positive for obliquity


# this would be the command to group by frequency
# note you must also make an equivalence for sigma -> abs(sigma) but also correct for te phase accordingly
# use package dplyr
# > OBLIQUITY_S <- OBLIQUITY %>% group_by(sigma) %>% summarise (N=sum(N))



PTERM1 <- list(  N = Dfseconde * outer(EPI$M,EPI$M), 
                       rate = outer(EPI$g, EPI$g, '-'), 
                       phase = outer(EPI$beta, EPI$beta, '-'))


PTERM2 <- list ( N = D1*IOM$N, rate = IOM$sigma, phase = IOM$delta)
PTERM3 <- list ( N = D2*IOM$N^2, rate = 2*IOM$sigma, phase = 2*IOM$delta)
PTERM4 <- list ( N = D2ik*outer(IOM$N, IOM$N), 
                       rate = outer(IOM$sigma, IOM$sigma, '+'), 
                       phase = outer(IOM$delta, IOM$delta, '+'))


PTERM5 <- list ( N = D2ikprime*outer(IOM$N, IOM$N), 
                       rate = outer(IOM$sigma, IOM$sigma, '-'), 
                       phase = outer(IOM$delta, IOM$delta, '-'))


for (i in seq(ng)) {
  for (k in seq(ng)) {
      if (k >= i)  {
        PTERM1$N[i, k] = NA;
        PTERM1$rate[i, k ] = NA;
        PTERM1$phase[i, k ] = NA;
      }
  }}

for (i in seq(nf)) {
  for (k in seq(nf)) {
      if (k >= i)  {
        PTERM4$N[i, k] = NA;
        PTERM4$rate[i, k ] = NA;
        PTERM4$phase[i, k ] = NA;
        PTERM5$N[i, k] = NA;
        PTERM5$rate[i, k ] = NA;
        PTERM5$phase[i, k ] = NA;
      }
  }}


PTERM1$N <- PTERM1$N[which(!is.na(PTERM1$N))]*1
PTERM1$rate <- PTERM1$rate[which(!is.na(PTERM1$rate))]
PTERM1$phase <- PTERM1$phase[which(!is.na(PTERM1$phase))]
PTERM4$N <- PTERM4$N[which(!is.na(PTERM4$N))]
PTERM4$rate <- PTERM4$rate[which(!is.na(PTERM4$rate))]
PTERM4$phase <- PTERM4$phase[which(!is.na(PTERM4$phase))]
PTERM5$N <- PTERM5$N[which(!is.na(PTERM5$N))]
PTERM5$rate <- PTERM5$rate[which(!is.na(PTERM5$rate))]
PTERM5$phase <- PTERM5$phase[which(!is.na(PTERM5$phase))]

PSI <- data.frame (
        N = c(PTERM1$N, PTERM2$N, PTERM3$N, PTERM4$N, PTERM5$N), 
        rate = c(PTERM1$rate, PTERM2$rate, PTERM3$rate, PTERM4$rate, PTERM5$rate), 
        phase = c(PTERM1$phase, PTERM2$phase, PTERM3$phase, PTERM4$phase, PTERM5$phase))



negatives <- which(PSI$rate < 0);

PSI$N[negatives] = -PSI$N[negatives]
PSI$rate[negatives] = -PSI$rate[negatives]
PSI$phase[negatives] = -PSI$phase[negatives]


require(dplyr)
PSI_ORDERED <- PSI %>% group_by(rate) %>% summarise (N=sum(N))

ORDER <- order(abs(PSI_ORDERED$N),  decreasing=TRUE)

PSI_ORDERED <- as.data.frame(PSI_ORDERED[ORDER, ])

plot(PSI_ORDERED$rate, abs(PSI_ORDERED$N), type='h')
plot(BRate, abs(BAmp), type='h')
