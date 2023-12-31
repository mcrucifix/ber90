dyn.load("p7505ff_SUBROUTINE.so")

pir =  pi/180
pirr = pir/3600

sectorad <- pi/(180*60.*60.)
degtorad <- pi/(180.)


  times <- seq(-1e6, 0, 1e3)
  require(palinsol)
  ORB <- sapply(times,ber90)
  ORB_LA <- sapply(times,la04)
  obl_ber90 <- ORB['eps',]*180/pi
  obl_la04 <- ORB_LA['eps',]*180/pi
  
  
  esinw_ber90 <- ORB['ecc',]*sin(-ORB['varpi',] )
  varpi_ber90 <- ORB['varpi',]
  

  esinw_la04 <- ORB_LA['ecc',]*sin(-ORB_LA['varpi',] )
  varpi_la04 <- ORB_LA['varpi',]
  
  e_ber90 <- ORB['ecc',]
 
compute_constants <- function(bea = 23.44579, prea = 50.273147, prega = -2514.27, ala = 54.9066, apoa = 17.3919 ) {
  result <- .Fortran("berger", 
                   bea = as.double ( bea ), 
                   prea = as.double ( prea ), 
                   prega = as.double ( prega ), 
                   ala = as.double ( ala ), 
                   apoa = as.double ( apoa ), 
                   prma = as.double(0),
                   pprma = as.double(0),
                   ha = as.double(0), 
                   tseta = as.double(0), 
                   aa = as.double(rep(0,80)), 
                   a = as.double(rep(0,80)), 
                   c = as.double(rep(0,80)), 
                   bb = as.double(rep(0,500)), 
                   b = as.double(rep(0,500)), 
                   d = as.double(rep(0,500)), 
                   bf = as.double(rep(0,65000)), 

  e_ber90 <- ORB['ecc',]
 
compute_constants <- function(bea = 23.44579, prea = 50.273147, prega = -2514.27, ala = 54.9066, apoa = 17.3919 ) {
  result <- .Fortran("berger", 
                   bea = as.double ( bea ), 
                   prea = as.double ( prea ), 
                   prega = as.double ( prega ), 
                   ala = as.double ( ala ), 
                   apoa = as.double ( apoa ), 
                   prma = as.double(0),
                   pprma = as.double(0),
                   ha = as.double(0), 
                   tseta = as.double(0), 
                   aa = as.double(rep(0,80)), 
                   a = as.double(rep(0,80)), 
                   c = as.double(rep(0,80)), 
                   bb = as.double(rep(0,500)), 
                   b = as.double(rep(0,500)), 
                   d = as.double(rep(0,500)), 
                   bf = as.double(rep(0,65000)), 
                   pf = as.double(rep(0,10000)), 
                   dpf = as.double(rep(0,10000)), 
                   sa = as.double(rep(0,10000)), 
                   ddr = as.double(rep(0,10000)))

  computed_constants = with(result, c(ktilde=pprma, eps=ha, k=prma, sigma=tseta) ) 
  return(list(computed_constants=computed_constants, full_results=result))

}


compute_solution <- function(result) {
  
  estar = result$ha
  psibar = result$prma * sectorad
  zeta   = result$tseta * degtorad
  
  ECC <- data.frame(
    V1 = seq(length(result$aa)),
    V2 = result$aa,
    V3 = result$a,
    V4 = result$c )
  
  EW <- data.frame(
    V1 = seq(length(result$aa)), 
    V2 = result$aa,
    V3 = result$a + result$prma,
    V4 = result$c + result$tseta)
  
  EW$V5 = 360*60*60 / EW$V3
  
  I2 <- data.frame( 
    V1 = seq(length(result$bb)), 
    V2 = result$bb,
    V3 = result$b,
    V4 = result$d)
  
  OBL <- data.frame(
    V1 = seq(10000),
    V2 = result$bf[seq(10000)] / pirr,     # 65000 rows
    V3 = result$pf / pirr,     # 10000 rows
    V4 = result$dpf / pirr,    # 10000 rows
    V5 = result$sa / pirr,     # 10000 rows
    V6 = result$ddr / pir)     # 10000 rows
  tables = list(I2=I2, OBL=OBL, ECC=ECC, EW=EW, estar=estar, psibar=psibar, zeta=zeta)
  return(tables)
} 


# pour retrouver exactement les parametres publies dans ber90
# et egalement encodes dans palinsol il faudrait utiliser: 
#        bea        prea       prega         ala        apoa
#     23.44580    50.27251 -2464.25259    54.88959    18.12131
# valeeurs originelles
#     23.44579    50.2731  -2514.27       54.9066     17.33919
# 

ber90_values<-c(bea = 23.44579, prea = 50.273147, prega = -2514.27, ala = 54.9066, apoa = 17.3919 )
optim_values<-c(bea=23.44579,prea=50.27251,prega= -2464.25259,ala=54.88959,apoa=18.12131)


result=compute_constants ( 
                    ber90_values['bea'], 
                    ber90_values['prea'], 
                    ber90_values['prega'], 
                    ber90_values['ala'], 
                    ber90_values['apoa'])$full_results

result=compute_constants ( 
                    optim_values['bea'], 
                    optim_values['prea'], 
                    optim_values['prega'], 
                    optim_values['ala'], 
                    optim_values['apoa'])$full_results





SOL = compute_solution(result)

test_solutions <- function(SOL, iplot=TRUE) {
  ECC=SOL$ECC
  I2=SOL$I2
  OBL=SOL$OBL
  EW=SOL$EW
  psibar=SOL$psibar
  zeta=SOL$zeta
  estar=SOL$estar

  indices <- order(OBL$V2)
  
  obl = sapply(times, function(t)  estar + sum(OBL$V2[indices]/3600*cos(OBL$V5[indices]*sectorad*t + OBL$V6[indices]*degtorad)))
  
  
 
  # e ber90 parfaitement reproduit
  
  # methode 1
  
  esinw = sapply(times, function(t)  sum(EW$V2*sin(EW$V3*sectorad*t + EW$V4*degtorad)))
  ecosw = sapply(times, function(t)  sum(EW$V2*cos(EW$V3*sectorad*t + EW$V4*degtorad)))
  
  varpi_method1 = Arg(ecosw + (0+1i)*esinw)
  
  # diff with ber90 of the order of 0.006 in 1 Myr
  
  # method 2 
  
  esinpi = sapply(times, function(t)  sum(ECC$V2*sin(ECC$V3*sectorad*t + ECC$V4*degtorad)))
  ecospi = sapply(times, function(t)  sum(ECC$V2*cos(ECC$V3*sectorad*t + ECC$V4*degtorad)))
  
  
  
#   obl = sapply(times, function(t)  estar + sum(OBL$V2/3600*cos(OBL$V5*sectorad*t + OBL$V6*degtorad)))
  
  
  psi = zeta +  psibar*times + sapply(times, function(t)  sum(OBL$V3*sectorad*sin(I2$V3*sectorad*t + I2$V6*degtorad)))

# zeta = 0.02793538
# psibar = 0.0002443016

# in ber90 palinsol

  psibar_ber90<- 50.41726176/60./60. * pi/180
  estar_ber90 <- 23.33340950
  zeta_ber90  <- 1.60075265 * pi/180.

e <- sqrt(esinpi^2+ecospi^2)


Pi <-atan(esinpi/ecospi)+pi*(ecospi<0)
varpi <- (Pi+psi+pi) %% (2*pi)

esinw_method2 <- e*sin(varpi)

# diff entre method1 et method2 de l'ordre de 0.002 (max)
# + difference de 180 degres sur varpi

if (iplot)
{
  plot(times, obl - obl_ber90, type='l')
  plot(obl, type='l', xlim=c(0,100))
  lines(obl_ber90, col='red')
  plot(times, e - e_ber90, type='l')
  plot(esinw, typ='l', xlim=c(0,100))
  lines(esinw_ber90, col='red')
  lines(esinw_la04, col='green', lwd=3)
  lines(times, esinw_ber90, typ='l', col='red')
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
}
out = sum((esinw_ber90 + esinw)^2)
print(out)
return(out)
}


# j'ai toujours des differences entre varpi (method1) et varpi (method2), de l'ordre de 0.02 et qui derivent
# attention aussi il y une differene 'pi' entre le deux
# avec la method1, je suis plus proche du ber90 officiel (calcule par palinsol avec la method2)
# a mon avis c'est parce que la methode 2 utilise beaucoup plus de termes
# il faudrait se forcer à n'utiliser que les 80 plus grands et voir ce que ça donne

# diff method2 en utilisant exactement le meme psibar que ber90: 0.002, ce qui est vraiment tr.es faible, et pas de dérive au cours du temps. 

# DONC: AVEC LE MEME PSIBAR que dans BER90, on obtient des differences de 0.002 sur esinw_method2 et 10e-4 sur obl.
# LA differences sur OBL pourrait tenir au nombre de termes retenus. 

print ( compute_constants()$computed_constants )

u=c(bea = 23.44579, prea = 50.273147, prega = -2514.27, ala = 54.9066, apoa = 17.3919)

to_optim <- function(u=c(bea = 23.44579, prea = 50.273147, prega = -2514.27, ala = 54.9066, apoa = 17.3919)){
#   u['ala']=54.9066
#   u['apoa']=17.3919
 out_t <- compute_constants(u['bea'], u['prea'], u['prega'], u['ala'], u['apoa'])
 result <- out_t$full_results
 SOL <- compute_solution(result)
 out <- out_t$computed_constants
 diff_target <- sum (with(as.list(out), c(ktilde - 50.417262, eps - 23.333410 , k - 50.390811, sigma - 1.600753 ) **2 ) ) 
 print(diff_target)
# return(test_solutions(SOL, iplot=FALSE))
 return(diff_target)
}

print(to_optim())

O = optim(u, to_optim)

