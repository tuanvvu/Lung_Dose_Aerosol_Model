#####ICRP MODELS FOR HYDROPHOBIC PARTICLES

###01 INPUT Diamter
# Input data
dth<-c(0,005,0.01,0.02,0.03,0.05,0.07,0.1,0.12,0.15,0.17,0.2,0.25,0.3,0.35,0.4,0.5,0.7,1,2,3,4,5,6,7,8,9,10,15,20)
# dae<-dth  ## Assuming the particles are sphere and desity is 1.
#dae<-c(0.0011,0.0021,0.0041,0.0101,0.019,0.038,0.090,0.168,0.312,0.5,0.7,1,2,3,5,7,10,15,20)

#### Convert dth(dB) to dae (da) diameter
dB<- 0.149    ### Inital Mobility diametter
da<-seq(0.01,20, 0.00005)  ### Aerodynamic diametter range from 0.01 um to 20 um
cA<-da
x<-cA
p<- 10 ### density measured by APS/SMPS algorithm.
X<-1 #### Shape factor
cB<- 1 + (0.0712/dB)*(2.514 + 0.800*exp(-0.55*(dB/0.0712))) ### Cunning factors for dB

y<-function(dB) {
for (i in 1:399801) {          
cA[i]<-sqrt((1 + (0.0712/da[i])*(2.514 + 0.800*exp(-0.55*(da[i]/0.0712)))))
x[i]<- abs(da[i] - sqrt(p/X)*dB*sqrt((1 + (0.0712/dB)*(2.514 + 0.800*exp(-0.55*(dB/0.0712)))))/cA[i])}
da[which(x == min(x), arr.ind = TRUE)] } ### from i=1 to length(da)###

MBD<- read.csv("Mobilitydiameter.csv",header=FALSE) ## dataframe of mobilitydiameter
dth<- as.numeric(as.vector(MBD[1,]))  ### Input Mobilitydiameter

AD<-MBD ### Aerodynamic diameter
for (j in 1:29) {                      
dB<-MBD[1,j]
AD[1,j]<-y(dB) }
dae<- as.numeric(as.vector(AD[1,]))   ##51<-leng(dth): number of sizebins##### Input Aerodynamicdiamter
# write.csv(dae,paste("C:/ICRP/","AerodynamicDiameterP2.csv",sep=""))

### 02. Subject Parameter

SFt=1 # for male adult
SFb=1 # for male adult
SFA=1 # for male adult
V= 300 #mL/s= cm3/s
Vn=300 # for male adult
VT= 750 #mL
VdBB= 49 
Vdbb= 47
VdET= 50
FRC= 3301
tB<- VdBB*(1+0.5*VT/FRC)/V
tb<- Vdbb*(1+0.5*VT/FRC)/V
tA<- (VT-VdET-(VdBB+Vdbb)*(1+VT/FRC))/V

###03. ICRP EQUATIONS FOR REGIONAL LUNG DEPOSITION
###03.1. Themordyamic deposition
# Empirical correction factor to allow for enhancement of thermodynamic deposition caused by nonlaminar bronchial airflow
Yth<- 1 + 100*exp(-(log10(100+10/(dth^0.9)))^2) 
# Slip correction
Lamda= 0.0712 # Mean free path (um)
Cd<-1 + Lamda/dth*(2.514+0.80*exp(-0.55*dth/Lamda))
# Diffusion coefficient
k=1.38*(10^-23) # Boltzman Constant, J/K= N.m/K
T=310  # Temperature (K), 37oC
miu=1.90*(10^-5) # Air viscosity, Pa.s= N.s/m2;  from Hindbook, miu=1.81E-05
#Diffusion coeffience 
D<- (k*(10^2))*T*Cd/(3*pi*(miu*(10^-4))*(dth*(10^-4))) #cm2/s
ath1= 18        ### a value for ET1 Inhalation during thermodynamic
ath2= 15.1     ### a value for ET2  Inhalation during thermodynamic
ath3= 22.02*(SFt^1.24)*Yth    ### a value for BB Inhalation during thermodynamic
ath4= -76.8 + 167*(SFb^0.65)  ### a value for bb Inhalation during thermodynamic
ath5= 170 + 103*(SFA^2.13)    ### a value for AL  
ath6= -76.8 + 167*(SFb^0.65)  ### a value for bb Exhalation
ath7=  22.02*(SFt^1.24)*Yth   ### a value for BB Enhalation
ath8= 15.1                    ### a value for ET2 Exhalation 
ath9= 18                      ### a value for ET1 Inhalation during thermodynmic
Rth1<-D*((Vn*SFt)^(-1/4))     ### R value for ET1 Inhalation during thermodynamic
Rth2<-D*((Vn*SFt)^(-1/4))   
Rth3<-D*tB
Rth4<-D*tb
Rth5<-D*tA
Rth6<-D*tb
Rth7<-D*tB
Rth8<-D*((Vn*SFt)^(-1/4)) 
Rth9<-D*((Vn*SFt)^(-1/4)) 
pth1=0.5                     ### p value for ET1 Inhalation during thermodynamic
pth2=0.538
pth3=0.6391
pth4=0.5676
pth5=0.6101
pth6=0.5676
pth7=0.6391
pth8=0.538
pth9=0.5
Q1=1                        ### Volumetric Fraction
Q2=1
Q3<- 1 -VdET/VT
Q4<- 1- (VdET + VdBB*(1+ VT/FRC))/VT
Q5<- 1 - (VdET+ VdBB*(1+ VT/FRC) + Vdbb*(1+VT/FRC))/ VT
Q6<- 1- (VdET + VdBB*(1+ VT/FRC))/VT
Q7<- 1 -VdET/VT
Q8=1
Q9=1
E2<-Q2/Q1                   ### The Quotient Volumetric Fraction
E3<-Q3/Q2
E4<-Q4/Q3
E5<-Q5/Q4
E6<-Q6/Q5
E7<-Q7/Q6
E8<-Q8/Q7
E9<-Q9/Q8
# Filter thermodynamic deposition efficiency
nth1<- 0.5*(1 - exp(-ath1*(Rth1^pth1)))
nth2<- 1 - exp(-ath2*(Rth2^pth2))
nth3<- 1 - exp(-ath3*(Rth3^pth3))
nth4<- 1 - exp(-ath4*(Rth4^pth4))
nth5<- 1 - exp(-ath5*(Rth5^pth5))
nth6<- 1 - exp(-ath6*(Rth6^pth6))
nth7<- 1 - exp(-ath7*(Rth7^pth7))
nth8<- 1 - exp(-ath8*(Rth8^pth8))
nth9<- 0.05*(1 - exp(-ath9*(Rth9^pth9)))

###03.2. Aerodynamic deposition
aae1=3.0*(10^-4)              ### a value for ET1 Inhalation during aerodynamic transport
aae2=5.5*(10^-5)
aae3=4.08*(10^-6)
aae4=0.1147
aae5=0.146*(SFA^0.98)
aae6=0.1147
aae7=4.08*(10^-6)
aae8=5.5*(10^-5)
aae9=3.0*(10^-4)
Rae1<-(dae^2)*Vn*(SFt^3)       ### R value for ET1 Inhalation during aerodynamic transport
Rae2<-(dae^2)*Vn*(SFt^3)
Rae3<-(dae^2)*V*(SFt^2.3)
Rae4<-(0.056 + (tb^1.5))*(dae^(tb^-0.25))
Rae5<-(dae^2)*tA
Rae6<-(0.056 + (tb^1.5))*(dae^(tb^-0.25))
Rae7<-(dae^2)*V*(SFt^2.3)
Rae8<-(dae^2)*Vn*(SFt^3)
Rae9<-(dae^2)*Vn*(SFt^3)
pae1=1                          ### p value for ET1 Inhalation during aerodynamic transport
pae2=1.17
pae3=1.152
pae4=1.173
pae5=0.6495
pae6=1.173
pae7=1.152
pae8=1.17
pae9=1
# Filter Aerodynamic deposition efficiency for regional lung 
nae1<- 0.5*(1 - 1/((aae1*(Rae1^pae1))+1))
nae2<- 1 - 1/((aae1*(Rae1^pae1))+1)
nae3<- 1 - exp(-aae3*(Rae3^pae3))
nae4<- 1 - exp(-aae4*(Rae4^pae4))
nae5<- 1 - exp(-aae5*(Rae5^pae5))
nae6<- 1 - exp(-aae6*(Rae6^pae6))
nae7<- 1 - exp(-aae7*(Rae7^pae7))
nae8<- 1 - 1/((aae9*(Rae8^pae9))+1)
nae9<- 0.5*(1 - 1/((aae9*(Rae9^pae9))+1))

# Deposition efficiency for the initial filter
nthi= 1-0.5*(1-(7.6*(10^-4)*(dth^2.8) +1)^(-1)) + 1*(10^-5)*exp(0.055*dth)
nth0=1-nthi ## For particle smaller than 1um, the nthi=0.5
naei= 1-0.5*(1-(7.6*(10^-4)*(dae^2.8) +1)^(-1)) + 1*(10^-5)*exp(0.055*dae) ### Ambient air velocity smaller than 1 m/s
nae0=1-naei
# For a calculation of windspeeds
# naei= 1-0.5*(1-(7.6*(10^-4)*(dae^2.8) +1)^(-1)) + 1*(10^-5)*(U^2.75)*exp(0.055*dae) #U: windspeed (m/s) from ICRP model
# For U< 4 m/s: naei= 0.5*(1+exp(-0.06*dae))  # From Hind's book 2002

## Deposition by themordynamic
DEth1<- nth1*(1-nth0)
DEth2<- DEth1*nth2*E2*(1/nth1 -1)
DEth3<- DEth2*nth3*E3*(1/nth2 -1)
DEth4<- DEth3*nth4*E4*(1/nth3 -1)
DEth5<- DEth4*nth5*E5*(1/nth4 -1)
DEth6<- DEth5*nth6*E6*(1/nth5 -1)
DEth7<- DEth6*nth7*E7*(1/nth6 -1)
DEth8<- DEth7*nth8*E8*(1/nth7 -1)
DEth9<- DEth8*nth9*E9*(1/nth8 -1)
DEthET<-DEth1 + DEth2 + DEth8 +DEth9
DEthTB<-DEth3 + DEth4 + DEth6 +DEth7
DEthAL<-DEth5
DEthtotal<-DEthET+DEthTB+DEthAL
## Deposition by Aerodynamic 
DEae1<- nae1*(1-nae0)
DEae2<- DEae1*nae2*E2*(1/nae1 -1)
DEae3<- DEae2*nae3*E3*(1/nae2 -1)
DEae4<- DEae3*nae4*E4*(1/nae3 -1)
DEae5<- DEae4*nae5*E5*(1/nae4 -1)
DEae6<- DEae5*nae6*E6*(1/nae5 -1)
DEae7<- DEae6*nae7*E7*(1/nae6 -1)
DEae8<- DEae7*nae8*E8*(1/nae7 -1)
DEae9<- DEae8*nae9*E9*(1/nae8 -1)
DEaeET<-DEae1 + DEae2 + DEae8 +DEae9
DEaeTB<-DEae3 + DEae4 + DEae6 +DEae7
DEaeAL<-DEae5
DEaetotal<-DEaeET+DEaeTB+DEaeAL


###03.3 COMBINING THERMODYNAMIC and AERODYNAMIC DEPOSITION
n1<- (nth1^2 + nae1^2)^(1/2)
n2<- (nth2^2 + nae2^2)^(1/2)
n3<- (nth3^2 + nae3^2)^(1/2)
n4<- (nth4^2 + nae4^2)^(1/2)
n5<- (nth5^2 + nae5^2)^(1/2)
n6<- (nth6^2 + nae6^2)^(1/2)
n7<- (nth7^2 + nae7^2)^(1/2)
n8<- (nth8^2 + nae8^2)^(1/2)
n9<- (nth9^2 + nae9^2)^(1/2)

DE1<- n1*(1-nae0) ## Assuming nth0=0 for small particles (Dp <100nm)
DE2<- DE1*n2*E2*(1/n1 -1)
DE3<- DE2*n3*E3*(1/n2 -1)
DE4<- DE3*n4*E4*(1/n3 -1)
DE5<- DE4*n5*E5*(1/n4 -1)
DE6<- DE5*n6*E6*(1/n5 -1)
DE7<- DE6*n7*E7*(1/n6 -1)
DE8<- DE7*n8*E8*(1/n7 -1)
DE9<- DE8*n9*E9*(1/n8 -1)

DEET<-DE1 + DE2 + DE8 +DE9
DETB<-DE3 + DE4 + DE6 +DE7
DEAL<-DE5
DEtotal<-DEET+DETB+DEAL

###04.OUTPUT THE RESULTS
write.csv(dae,paste(workingDirectory,"daeP10mansit.csv",sep=""))
write.csv(DEET,paste(workingDirectory,"DEETPP10mansit.csv",sep=""))
write.csv(DETB,paste(workingDirectory,"DETBPP10mansit.csv",sep=""))
write.csv(DEAL,paste(workingDirectory,"DEAPLP10mansit.csv",sep=""))
write.csv(DEtotal,paste(workingDirectory,"DETotalP10mansit.csv",sep=""))
