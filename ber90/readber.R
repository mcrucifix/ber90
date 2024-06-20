filename = 'ber907505_1'


sectorad <- pi/(180*60.*60.)
degtorad <- pi/(180.)

L1  = read.table(filename, skip=3, nrow=1)

Eps0 = L1$V1
k0 = L1$V3
k1 = L1$V4

Eps0

L2  = read.table(filename, skip=4, nrow=1)

n = L2$V1
m = L2$V2
no = L2$V3
estar = L2$V5
psibar = L2$V6 * sectorad
zeta = L2$V8* degtorad

ECC = read.table(filename, skip=5, nrow=80);
I2  = read.table(filename, skip=85, nrow=80);
OBL = read.table(filename, skip=165, nrow=6480);
EW  = read.table(filename, skip=9805, nrow=212);

# columns for ECC are (Ampritude, mean rate, phase)
# columns for EW are amplitude, mean rate, hpase, period
# columns for OBL are 

# obliquity :  numeber, amplitude is in seconds; ?, ?, rate in ''/year, phase in degrees, 'bp', '..'

# OBLIQUITY

times <- seq(-1e6, 0, 1e3)

indices <- order(OBL$V2)

obl = sapply(times, function(t)  estar + sum(OBL$V2[indices]/3600*cos(OBL$V5[indices]*sectorad*t + OBL$V6[indices]*degtorad)))
plot(obl)


require(palinsol)
ORB <- sapply(times,ber90)
obl_ber90 <- ORB['eps',]*180/pi
plot(times, obl - obl_ber90, type='l')
plot(obl, type='l', xlim=c(0,100))
lines(obl_ber90, col='red')


esinw_ber90 <- ORB['ecc',]*sin(-ORB['varpi',] )
varpi_ber90 <- ORB['varpi',]

e_ber90 <- ORB['ecc',]

# e ber90 parfaitement reproduit

# methode 1

esinw = sapply(times, function(t)  sum(EW$V2*sin(EW$V3*sectorad*t + EW$V4*degtorad)))
ecosw = sapply(times, function(t)  sum(EW$V2*cos(EW$V3*sectorad*t + EW$V4*degtorad)))

varpi_method1 = Arg(ecosw + (0+1i)*esinw)

# diff with ber90 of the order of 0.006 in 1 Myr

# method 2 

esinpi = sapply(times, function(t)  sum(ECC$V2*sin(ECC$V3*sectorad*t + ECC$V4*degtorad)))
ecospi = sapply(times, function(t)  sum(ECC$V2*cos(ECC$V3*sectorad*t + ECC$V4*degtorad)))



obl = sapply(times, function(t)  estar + sum(OBL$V2/3600*cos(OBL$V5*sectorad*t + OBL$V6*degtorad)))


psi = zeta +  psibar*times + sapply(times, function(t)  sum(OBL$V3*sectorad*sin(I2$V3*sectorad*t + I2$V6*degtorad)))

# zeta = 0.02793538
# psibar = 0.0002443016

# in ber90 palinsol

psibar_ber90<- 50.41726176/60./60. * pi/180
estar_ber90 <- 23.33340950
zeta_ber90  <- 1.60075265 * pi/180.




e <- sqrt(esinpi^2+ecospi^2)

plot(times, e - e_ber90, type='l')

Pi <-atan(esinpi/ecospi)+pi*(ecospi<0)
varpi <- (Pi+psi+pi) %% (2*pi)

esinw_method2 <- e*sin(varpi)

# diff entre method1 et method2 de l'ordre de 0.002 (max)
# + difference de 180 degres sur varpi

plot(times, esinw, typ='l')
plot(times, esinw_ber90 - esinw, col='red')
# differnces de 0.006, grandissantes au cours de temps

plot(times, esinw_ber90 + esinw_method2, col='red')
# differnces de 0.005, avec u pic .a 0.010 à -1Myr
# plus difference de 180 degres sur la definition de pi

plot(Pi)
lines(Arg(ecospi + esinpi*(0+1i)))
lines(Arg(ecospi + esinpi*(0+1i)), col='red')

plot(varpi, type='l', xlim=c(900,1000))
lines(varpi_ber90, col='red')


# j'ai toujours des differences entre varpi (method1) et varpi (method2), de l'ordre de 0.02 et qui derivent
# attention aussi il y une differene 'pi' entre le deux
# avec la method1, je suis plus proche du ber90 officiel (calcule par palinsol avec la method2)
# a mon avis c'est parce que la methode 2 utilise beaucoup plus de termes
# il faudrait se forcer à n'utiliser que les 80 plus grands et voir ce que ça donne

# diff method2 en utilisant exactement le meme psibar que ber90: 0.002, ce qui est vraiment tr.es faible, et pas de dérive au cours du temps. 

# DONC: AVEC LE MEME PSIBAR que dans BER90, on obtient des differences de 0.002 sur esinw_method2 et 10e-4 sur obl.
# LA differences sur OBL pourrait tenir au nombre de termes retenus. 


