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
I2  = read.table(filename, skip=75, nrow=80);
OBL = read.table(filename, skip=165, nrow=6480);
EW  = read.table(filename, skip=9805, nrow=212);

# columns for ECC are (Ampritude, mean rate, phase)
# columns for EW are amplitude, mean rate, hpase, period
# columns for OBL are 

# obliquity :  numeber, amplitude is in seconds; ?, ?, rate in ''/year, phase in degrees, 'bp', '..'

# OBLIQUITY

times <- seq(-1e6, 0, 1e3)



obl = sapply(times, function(t)  estar + sum(OBL$V2/3600*cos(OBL$V5*sectorad*t + OBL$V6*degtorad)))
plot(obl)


require(palinsol)
ORB <- sapply(times,ber90)
obl_ber90 <- ORB['eps',]*180/pi

plot(obl - obl_ber90, type='l', xlim=c(0,10))
plot(obl, type='l', xlim=c(0,100))
lines(obl_ber90, col='red')


esinw_ber90 <- ORB['ecc',]*sin(-ORB['varpi',] )
varpi_ber90 <- ORB['varpi',]

e_ber90 <- ORB['ecc',]

# methode 1

esinw = sapply(times, function(t)  sum(EW$V2*sin(EW$V3*sectorad*t + EW$V4*degtorad)))
ecosw = sapply(times, function(t)  sum(EW$V2*cos(EW$V3*sectorad*t + EW$V4*degtorad)))

varpi_method1 = Arg(ecosw + (0+1i)*esinw)

# method 2 

esinpi = sapply(times, function(t)  sum(ECC$V2*sin(ECC$V3*sectorad*t + ECC$V4*degtorad)))
ecospi = sapply(times, function(t)  sum(ECC$V2*cos(ECC$V3*sectorad*t + ECC$V4*degtorad)))



obl = sapply(times, function(t)  estar + sum(OBL$V2/3600*cos(OBL$V5*sectorad*t + OBL$V6*degtorad)))

psi = zeta +  psibar*times + sapply(times, function(t)  sum(OBL$V3*sectorad*sin(I2$V3*sectorad*t + I2$V6*degtorad)))


e <- sqrt(esinpi^2+ecospi^2)
Pi <-atan(esinpi/ecospi)+pi*(ecospi<0)
varpi <- (Pi+psi+pi) %% (2*pi)


plot(times, esinw, typ='l', xlim=c(-5e3,0), ylim=c(0.00, 0.02))
lines(times, esinw_ber90, col='red')

plot(esinw - esinw_ber90, type='l')

plot(times, varpi - varpi_ber90, type='l', ylim=c(-0.15, 0.15))

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


