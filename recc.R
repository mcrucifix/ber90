
require(gtools)
EPI <- read.table('epi.dat')
colnames(EPI) <- c('index','omega','M','phi')


## Il y asans doute une erreur dans la page 54 de Berger et Loutre, mais ou ???

## Revoir la these de Berger ? 

print('init tables')

EPI$omega <- EPI$omega/(60*60*360)*2*pi
EPI$M  <- EPI$M / 1.e8
EPI$phi <- EPI$phi * 2*pi / 360. 

combine <- function(a,b, om1, om2, ph1, ph2) {
  A <- outer(a,b, "*")
  O <- outer(om1,om2, "+")
  P <- outer(ph1,ph2, "+")
  A <- A[upper.tri(A)]
  O <- O[upper.tri(O)]
  P <- P[upper.tri(P)]
  return(data.frame(A=A,O=O,P=P))
}

# order  by decreasing frequency

O1 <- order(EPI$omega, decreasing=TRUE)

EPI <- list(omega=EPI$omega[O1], M=EPI$M[O1], phi=EPI$phi[O1])

m <- sqrt(sum(EPI$M^2))
EPI$M = EPI$M / m


n=0
IJB = matrix(0,3160,2)
for (i in seq(79)) for (j in seq(i+1,80))  {n=n+1 ; IJB[n,]=c(i,j)}

# contsuct the gamma beta


B <- combine(EPI$M, EPI$M, EPI$omega, -EPI$omega,  EPI$phi, -EPI$phi)



bk <- B$A
gk <- B$O
pk <- B$P

# bk <- EPI$M[IJB[,1]] * EPI$M[IJB[,2]]
# gk <- EPI$omega[IJB[,1]] - EPI$omega[IJB[,2]]
# pk <- EPI$phi[IJB[,1]] - EPI$phi[IJB[,2]]

ek <- bk + 0.375*bk*bk*bk 
sumbl2 <- rep(sum(bk^2), length(bk))
sumbl2 <- sumbl2 - bk*bk

ek <- ek + 0.75 * bk*sumbl2

print('group 1 : termes d ordre 1')

group0 <- data.frame(A=1-0.25*sum(bk^2), O=0, P=0)

 group1 <- data.frame(A=ek, O=gk, P=pk)
#  group2 <- -0.25 * data.frame(A=bk^2, O=2*gk, P=2*pk)
#  group5 <- 0.125 * data.frame(A=bk^3, O=3*gk, P=3*pk)

 odecreasing <- order(bk, decreasing=TRUE)
 bk <- bk[odecreasing]
 gk <- gk[odecreasing]
 pk <- pk[odecreasing]

# print('groupes 3 a 9 : termes d ordre 2')

# n=0
# IJ = matrix(0, 124750,2)
# for (i in seq(499)) for (j in seq(i+1,500))  {n=n+1 ; IJ[n,]=c(i,j)}

# bk3 = -0.5 * (bk[IJ[,1]] * bk[IJ[,2]] )
# gk3a = (gk[IJ[,1]] + gk[IJ[,2]] )
# gk3b = (gk[IJ[,1]] - gk[IJ[,2]] )

# pk3a = (pk[IJ[,1]] + pk[IJ[,2]] )
# pk3b = (pk[IJ[,1]] - pk[IJ[,2]] )

# group3 =  data.frame(A=c(bk3,bk3), O=c(gk3a,gk3b), P=c(pk3a,pk3b))



# bk6 = 0.375 * (bk[IJ[,1]] * bk[IJ[,2]]^2 )
# bk7 = 0.375 * (bk[IJ[,1]]^2 * bk[IJ[,2]] )

# gk6a = (gk[IJ[,1]] + 2*gk[IJ[,2]] )
# gk6b = (-gk[IJ[,1]] + 2*gk[IJ[,2]] )

# pk6a = (pk[IJ[,1]] + 2*pk[IJ[,2]] )
# pk6b = (-pk[IJ[,1]] + 2*pk[IJ[,2]] )

# gk7a = (2*gk[IJ[,1]] + gk[IJ[,2]] )
# gk7b = (2*gk[IJ[,1]] - gk[IJ[,2]] )

# pk7a = (2*pk[IJ[,1]] + pk[IJ[,2]] )
# pk7b = (2*pk[IJ[,1]] - pk[IJ[,2]] )


# group6 =  data.frame(A=c(bk6,bk6), O=c(gk6a,gk6b), P=c(pk6a,pk6b))
# group7 =  data.frame(A=c(bk7,bk7), O=c(gk7a,gk7b), P=c(pk7a,pk7b))

group2 <- -0.25 * data.frame(A=bk^2, O=2*gk, P=2*pk)
group3 <- -0.5 * combine(bk, bk, gk, gk, pk, pk)
group4 <- -0.5 * combine(bk,bk, gk, -gk, pk, -pk)
group5 <- 0.125 * data.frame(A=bk^3, O=3*gk, P=3*pk)

group6 <- 0.375 * combine(bk^2, bk, 2*gk, gk, 2*pk, pk)
group7 <- 0.375 * combine(bk, bk^2, gk, 2*gk, pk, 2*pk)
group8 <- 0.375 * combine(bk^2, bk, 2*gk, -gk, 2*pk, -pk)
group9 <- 0.375 * combine(bk, bk^2, -gk, 2*gk, -pk, 2*pk)




IJK = gtools::combinations(50,3)

bklm =  0.75 * as.numeric(apply(IJK, 1, function(V)  {bk[V[1]] *  bk[V[2]] * bk[V[3]]}))
gklm1 = as.numeric(apply(IJK, 1, function(V) {gk[V[1]] +  gk[V[2]] + gk[V[3]]}))
gklm2 = as.numeric(apply(IJK, 1, function(V) {gk[V[1]] +  gk[V[2]] - gk[V[3]]}))
gklm3 = as.numeric(apply(IJK, 1, function(V) {gk[V[1]] -  gk[V[2]] + gk[V[3]]}))
gklm4 = as.numeric(apply(IJK, 1, function(V) {gk[V[1]] -  gk[V[2]] - gk[V[3]]}))
pklm1 = as.numeric(apply(IJK, 1, function(V) {pk[V[1]] +  pk[V[2]] + pk[V[3]]}))
pklm2 = as.numeric(apply(IJK, 1, function(V) {pk[V[1]] +  pk[V[2]] - pk[V[3]]}))
pklm3 = as.numeric(apply(IJK, 1, function(V) {pk[V[1]] -  pk[V[2]] + pk[V[3]]}))
pklm4 = as.numeric(apply(IJK, 1, function(V) {pk[V[1]] -  pk[V[2]] - pk[V[3]]}))

group10 =  data.frame(A=-c(bklm,bklm,bklm,bklm), 
                      O=c(gklm1,gklm2,gklm3,gklm4),
                      P=c(pklm1,pklm2,pklm3,pklm4))

print('developpement final')

developpement <- rbind(group0, group1, group2, group3, group4,  group5, group6, group7, group8, group9,  group10)

print(length(developpement$A))

# developpement <- group10
# ici on prend un petit risque, on ne va garder que les termes dont l'amplitude est > 5.e-7, soit environ 1000 termes. 

print('ecremage')


developpement <- developpement[which(abs(developpement$A) > 1.e-4),]

print(length(developpement$A))

print ('chasse aux doublons')

developpement_1 <- data.frame(A=developpement$A,  O=developpement$O)
developpement_2 <- data.frame(P=developpement$P,  O=developpement$O)

AFF1 <- aggregate(.~O, data=developpement_1, FUN=sum)
AFF2 <- aggregate(.~O, data=developpement_2, FUN=function(x) x[1])

AFF <- cbind(AFF1, P=AFF2$P)

print ('liste finale:')

print(length(developpement$A))


print('generation series temporelles')

times <- seq(0,2e3,2)*1e3
require(palinsol)

AFF$A = AFF$A*m
developpement$A = developpement$A*m


ecc_reconstruct_non_aggregated <- sapply(times, function(t) sum(developpement$A * cos(developpement$O*t+developpement$P)))

# p_reconstruct <- matrix(0, ncol=length(times), nrow=8)

# p_reconstruct[1,] <- sapply(times, function(t) sum(m*group0$A * cos(group0$O*t+group0$P)))
# p_reconstruct[2,] <- sapply(times, function(t) sum(m*group1$A * cos(group1$O*t+group1$P)))
# p_reconstruct[3,] <- sapply(times, function(t) sum(m*group2$A * cos(group2$O*t+group2$P)))
# p_reconstruct[4,] <- sapply(times, function(t) sum(m*group3$A * cos(group3$O*t+group3$P)))
# p_reconstruct[5,] <- sapply(times, function(t) sum(m*group5$A * cos(group5$O*t+group5$P)))
# p_reconstruct[6,] <- sapply(times, function(t) sum(m*group6$A * cos(group6$O*t+group6$P)))
# p_reconstruct[7,] <- sapply(times, function(t) sum(m*group7$A * cos(group7$O*t+group7$P)))
# p_reconstruct[8,] <- sapply(times, function(t) sum(m*group10$A * cos(group10$O*t+group10$P)))

esinpi <- sapply(times, function(t) sum(EPI$M * sin(EPI$omega*t+EPI$phi)))
ecospi <- sapply(times, function(t) sum(EPI$M * cos(EPI$omega*t+EPI$phi)))

e2 = sqrt(esinpi*esinpi+ecospi*ecospi) * m

##

print('ici')
print(length(AFF$A))

print('reconstruction...')
ecc_reconstruct <- sapply(times, function(t) sum(AFF$A * cos(AFF$O*t+AFF$P)))
print('palinsol.')
ecc_ber90 <- sapply(times, ber90)['ecc',]





plot(ecc_reconstruct, type='l')
# lines(ecc_reconstruct_non_aggregated, type='l', col='blue')
lines(ecc_ber90-0.02, col='red')

plot(ecc_reconstruct - ecc_ber90, type='l')

print('true variances')
# print(var(ecc_reconstruct_non_aggregated))
print(var(ecc_reconstruct))
print(var(ecc_ber90))
print(var(ecc_ber90 - ecc_reconstruct))

# Order1 <- order(AFF$O, decreasing=FALSE)
# Order2 <- order(AFF$A, decreasing=TRUE)
# 
# plot(2*pi/AFF$O[Order1], AFF$A[Order1], type='l', xlim=c(0,3e6))
# AFF$period = 2*pi/AFF$O
# 
# AFF[Order2[seq(10)],]
# 


