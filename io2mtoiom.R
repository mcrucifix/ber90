I2OM <- read.table('iom.dat')
colnames(I2OM) <- c('index','sigmaprime','Nprime','deltaprime')


I2OM <- I2OM[-1,]  # suppress constant term 

# need to think more about this. 

# transform into radians / year and radians

I2OM$sigmaprime <- I2OM$sigmaprime/(60*60*360)*2*pi
I2OM$Nprime  <- I2OM$Nprime / 1.e8
I2OM$deltaprime <- I2OM$deltaprime * 2*pi / 360. 


n <- length(I2OM$Nprime)

# I will privilege here transparent writing, even if less efficient. Fortran compiler would be better

# line 1

TERM1 <- list()
TERM1$sigma <- I2OM$sigmaprime
TERM1$N <- I2OM$Nprime
TERM1$delta <- I2OM$deltaprime

for (i in seq(n)) TERM1$N[i] <- TERM1$N[i] * ( 2 - I2OM$Nprime[i]^2 - 2 * sum(I2OM$Nprime[-i]^2) ) 


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

tol <- 1.e-7

nkeep <- which (IOM$N > tol)

IOM <- list( N = IOM$N[nkeep], 
             sigma = IOM$sigma[nkeep], 
             delta = IOM$delta[nkeep]) 



#################################

### ATTENTION: ARNAQUE !!

### EN FAIT CLAIREMENT DANS IOM C'est les angles pleins, pas les demis
### QUI SONT FOURNIS. DONC CE QUI EST AU DESSUS EST SUPERFLU

#################################

IOM <- I2OM

# no duplicates. 

#### development of obliquity

epsilonbar = 23.4 * pi / 180
psibar     = 50.273147 * 2 * pi / 60 / 360 / 60
zeta       = 0.02793538

# thesis Berger p. 110 = Sharaf Budnikova

fi <-  psibar + IOM$sigma
ci <-  psibar / fi 
deltaprime <- IOM$delta + zeta 

tane = tan(epsilonbar)
cote = 1./tane

nf <- length(ci)

C1f <- -ci
C2f <- -0.5 * outer (seq(nf), seq(nf), function(i,k) {
                       psibar/(fi[i]+fi[k])*((ci[i]^2 + ci[k]^2 - 2) * tane +
                              (ci[i]+ci[k]) * cote)  })


C2fii = -.25 *(ci * (ci^2 - 1) * tane + ci^2 * cote )

C2fprime <- -0.5 * outer (seq(nf), seq(nf), function(i,k) {
                       psibar/(fi[i]-fi[k])*((ci[i]^2 - ci[k]^2 + 2*ci[k] - 2*ci[i]) * tane +
                              (ci[i]-ci[k]) * cote)  })



## ATTENTION CECI DOIT ETRE CORRIGE
## C'EST EN FONCTION DE CIF

## mais je dois revoir la these pasge 113
## et convertir les cif en ci. 

## pour ca il faut calculer les Di, qui ne sont pas encore calcules !!!! 

### DONC IL Y A ENCORE DU TAF ICI
### EQUATTIONS 63 DE BERGER

for (i in seq(nf)) { C2fprime[i,i] = 0 }

C1 = C1f + 1


ETERM1 <- list(  N = C1 * IOM$N, sigma = fi, delta = deltaprime ) 

ETERM2 <- list ( N = C2fii * IOM$N^2 , sigma = 2*fi, delta = 2 * deltaprime ) 

ETERM3 <- list ( N = C2f * outer(IOM$N, IOM$N), sigma = outer(fi, fi, '+'), delta = outer(deltaprime, deltaprime, '+'))

for (i in seq(length(IOM$N))) {
  for (k in seq(length(IOM$N))) {
     if ( i >= k ) {
       ETERM3$N[i,k] = NA; 
       ETERM3$sigma[i,k] = NA;
       ETERM3$delta[i,k] = NA}}}


ETERM3$N <- ETERM3$N[which(!is.na(ETERM3$N))]
ETERM3$sigma <- ETERM3$sigma[which(!is.na(ETERM3$sigma))]
ETERM3$delta <- ETERM3$delta[which(!is.na(ETERM3$delta))]


ETERM4 <- list ( N = C2fprime * outer(IOM$N, IOM$N), sigma = outer(fi, fi, '-'), delta = outer(deltaprime, deltaprime, '-'))

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

nkeep = 500
OR <- order(abs(OBLIQUITY$N), decreasing=TRUE)[seq(nkeep)]
OBLIQUITY <- with(OBLIQUITY, list(N = N[OR], sigma=sigma[OR], delta=delta[OR]))

# with this, there seems be no duplicates (but maybe there will be once we have recitfied the N's)
# however we have a problem becase the sigma0 generates periods = psibar

which(duplicated(OBLIQUITY$sigma))

2*pi/ OBLIQUITY$sigma[seq(100)]

require(palinsol)
data(BER90)
par(mfrow=c(1,2))

plot(OBLIQUITY$sigma[seq(100)]/(2*pi), OBLIQUITY$N[seq(100)]/(2*pi), type='n')
lines(ETERM1$sigma/(2*pi), ETERM1$N/(2*pi), col='red', type='h')
#lines(ETERM2$sigma/(2*pi), ETERM2$N/(2*pi), col='blue' , type='h')
#lines(ETERM3$sigma/(2*pi), ETERM3$N/(2*pi), col='gray', type='h')
#lines(ETERM4$sigma/(2*pi), ETERM4$N/(2*pi), col='green', type='h')
with(BER90$Table1,  plot(Rate/360/3600 , Amp/60/3600, type='h', col='black'))

#plot(ETERM1$sigma/(2*pi), C1f, type='h')

# attention at this point we haven't tried to reorganise and chase doublons. 
# there will be many of those, but one step at a time. 





#### DEVELOPPEMENT FOR PRECESSION

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




Dif <- ci * (cote + (ci - 1) * tane)

Diif <- 0.25 * ci * (ci^2 + ci - 1) * 0.125 * ci^2*(ci-1)^2*(tane^2) + .5 * ci^2 * cote^2 

Dikf  <- outer(seq(nf), seq(nf) function(i,j){
  psibar/(fi[i]+fi[k])*(
                        0.5*(ci[i]^2 + ci[k]^2 + ci[i] + ci[k] - ci[i]*ci[k] - 2) 
                      + C2f[i,k]*tane
                      - 0.5*(ci[i]*(ci[i]-1) + ci[k] * (ci[k] - 1)) * tane^2 + (ci[i] + ci[k]) * cote^2  ) } )



# must define here P0, ell but that's another part of the code. 

Dfseconde = 3 * psibar * P0 / ell * outer(seq(ng), seq(ng), function(i,k) { 1./(EPI$g[i]-EPI$g[k])})


# psi still need to be coded. eq. 64 of his thesis p. 112


