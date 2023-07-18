# clear all variables in the workspace
rm(list=ls())

#clear the console
cat("\014")

library(latex2exp)

######################################
#        Egg  ----> Pupae            #
######################################
GeL <- function(Temp){
  
  a  <- 1.5909 *(10^(-4))
  TL <- 2.0919
  TM <- 32.088
  m  <- 4.0
  
  if( (Temp > TL) & (Temp < TM) ){
    
    out  <- a*Temp*(Temp - TL)*((TM - Temp)^(1/m))  
    
  }else{
    
    out  <- 0.001
    
  }
  
  return(out)
  
}


######################################
#        Pupae  ----> Adult          #
######################################
GP <- function(Temp){
  
  a  <- 2.3699 *(10^(-4))
  TL <- 4.0
  TM <- 33.164
  m  <- 4.0
  
  if( (Temp > TL) & (Temp < TM) ){
    
    out  <- a*Temp*(Temp - TL)*((TM - Temp)^(1/m))  
    
  }else{
    
    out  <- 0.001
    
  }
  
  return(out)
  
}

######################################
#           Adult Survival           #
######################################
GA <- function(Temp){
  
  a  <- 6.8417 *(10^(-4))
  TL <- -3.0
  TM <- 30.034
  m  <- 2.50
  
  if( (Temp > TL) & (Temp < TM) ){
    
    out  <- a*Temp*(Temp - TL)*((TM - Temp)^(1/m))  
    
  }else{
    
    out  <- 0.001
    
  }
  
  return(out)
  
}

######################################
#           Mortality                #
######################################
Mortality <- function(Temp){
  
   # a1 <-  0.000002
   # b1 <- -0.000161
   # c1 <-  0.004191
   # d1 <- -0.046696
   # e1 <-  0.20067
   # 
   # out  <- (a1*(Temp^4)) + (b1*(Temp^3)) + (c1*(Temp^2)) + (d1*(Temp)) + e1
  
  k    <- 1.0
  TMAX <- 23.4265085
  rho  <- -5.5455493

  AA <-  (TMAX - Temp)/rho

  out <- 1.0 - (k*exp(1.0 + AA  - exp(AA) ) )
  
  if(out > 1.0) {
    
    out <- 1.0
    
  }

  return(out)
  
}


######################################
#           Fertility                #
######################################
beta <- function(Temp){
  
  #Average eggs /  females / days at optimal temperature  
  N <- 20
  
   alpha  <- 659.06
   gamma  <- 88.53
   lambda <- 52.32
   delta  <- 6.06
   tau    <- 22.87
   Tlow   <- 5.0
   Tmax   <- 30.0

  a <- (alpha*(gamma + 1))/(pi*(lambda^((2*gamma)+2)))

  b <- ( (lambda^2) - ((Temp - tau)^2) - (delta^2) )^gamma

  if( (Temp > Tlow) & (Temp < Tmax) ){

    out  <- N*a*b

  }else{

    out  <- 0.0

  }

  return(out)


  
#   #Female fecundity rate at diffeent temperature
#   a <-  20.875
#   b <-  8.125
#   
# 
#  # if( (Temp > Tlow) & (Temp < Tmax) ){
# 
#     out  <-  N*( 1.0 -  ( ((Temp - a)/b)^2 ) ) 
#     
# #    out  <- max(0.0,out)
# 
# # }else{
# 
#     out  <- max(0.0,out)
# 
# #  }
# 
#   return(out)
  
}

#################
#   Risk index  #
#################
RiskIndex <-function(Temp){
  
  ######################################
  #           Sex ratio                #
  ######################################
  SR <- 0.5
  
  M  <- Mortality(Temp)
  
  MP <- M
  
  MA <- M
  
  #print(SR)
  
  Numerator   <- beta(Temp)*SR*GP(Temp)*(GeL(Temp)^4)*(1-GA(Temp)-MA)
  
  Denominator <- ((GeL(Temp) + M)^4)*(GP(Temp) + MP)*(GA(Temp) + MA)
  
  out <-  Numerator/Denominator
  
  return(out)
  
}

Temp <-  seq(10,40,by=0.1)

Threshold <- rep(1,length(Temp))

RI <- rep(0,length(Temp))


for (i in 1:length(Temp)) {
  
  RI[i] <- RiskIndex(Temp[i])
  
}

RI

plot(Temp, 
     RI, 
     #lty=1,
     type="l",
     lwd = 4,
     xlab= TeX('$Temperature (C)$'), 
     ylab= "Risk Index",
     #xlim =c(0,50),
     ylim =c(0,3),
     col="blue",
     #pch = 16, 
     #cex=2,
     cex.lab=1.5, 
     cex.axis=1.5, 
     frame.plot = FALSE,
     axes=FALSE,
     font.lab=2)

axis(side = 1, lwd = 3 ,cex.lab=1.5, cex.axis =2, font.lab=2)
axis(side = 2, lwd = 3 ,cex.lab=1.5, cex.axis =2, font.lab=2)

lines(Temp,
      Threshold, 
      type="l", 
      lty=3, 
      lwd=4, 
      col="red")

legend("topright",
       legend=c("Risk Index ","Stability Threshold"),
       col=c("blue","red"),
       lty=c(1,3),
       lwd=4,
       ncol=1, cex = 1.5,text.font=2)

# Create a title with a red, bold/italic font
title(main="Risk Index for Drosophila suzukii", col.main="black", font.main=2, cex.main = 2)

