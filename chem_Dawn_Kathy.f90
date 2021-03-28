SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
IMPLICIT NONE
INTEGER NDIM, IJAC, ICP(*)
DOUBLE PRECISION U(NDIM), PAR(*), F(NDIM), DFDU(*), DFDP(*)
DOUBLE PRECISION kon1, koff1, kon2, koff2, kfs, L, n, B, I_rho, L_rho, &
delta_rho, L_R, I_R, delta_R, alpha_R, delta_P, I_K, k_X, k_G, k_C, PIX, &
Paxtot, m, alpha_PAK, Rho_Square, Rac_Square, Pax_Square, PAKtot, Pax, FAK, &
PaxFAK, Paxs, PaxsFAK, GIT, PaxGIT, PaxsGIT, Raci, Rac, Rhoi, Rho, cnsrv_1, &
cnsrv_2, cnsrv_3, cnsrv_4, cnsrv_5, R, RacRatio, RhoRatio, PaxRatio, K_is, &
K0, K, I_Ks, P_i, Rbar, Q_R, Q_rho, serine_kinase, Q_ps, J_gamma_2_1, &
J_gamma_1_4, J_gamma_2_4, J_gamma_1_6, J_gamma_2_6, J_gamma_1_10, &
J_gamma_2_10, f_FAK, f_PaxFAK, f_Paxs, f_GIT, f_PaxGIT, f_Raci, f_Rhoi

PaxFAK=U(1)
Paxs=U(2)
PaxsFAK=U(3)
PaxGIT=U(4)
PaxsGIT=U(5)
Rac=U(6)
Rho=U(7)

kon1=PAR(1)
koff1=PAR(2)
kon2=PAR(3)
koff2=PAR(4)
kfs=PAR(5)
L=PAR(6)
n=PAR(7)
B=PAR(8)
I_rho=PAR(9)
L_rho=PAR(10)
delta_rho=PAR(12)
L_R=PAR(13)
I_R=PAR(14)
delta_R=PAR(15)
alpha_R=PAR(16)
delta_P=PAR(17)
I_K=PAR(18)
k_X=PAR(19)
k_G=PAR(20)
k_C=PAR(21)
PIX=PAR(22)
Paxtot=PAR(23)
m=PAR(24)
alpha_PAK=PAR(25)
Rho_Square=PAR(26)
Rac_Square=PAR(27)
Pax_Square=PAR(28)
PAKtot=PAR(29)

cnsrv_1=1
cnsrv_2=50
cnsrv_3=2
cnsrv_4=1
cnsrv_5=Rac_Square
Pax=(Paxs*Pax_Square + PaxFAK*Pax_Square + PaxGIT*Pax_Square + PaxsFAK*&
Pax_Square + PaxsGIT*Pax_Square - Pax_Square*1)/(K_is*PAKtot*PIX*PaxGIT*k_C*&
k_G*k_X - Pax_Square + K_is*PAKtot*PIX*PaxsGIT*k_C*k_G*k_X - K_is*PAKtot*PIX*&
2*k_C*k_G*k_X + K_is*PAKtot*PIX*PaxGIT*R*alpha_R*k_C*k_G*k_X + K_is*PAKtot*&
PIX*PaxsGIT*R*alpha_R*k_C*k_G*k_X - K_is*PAKtot*PIX*R*alpha_R*2*k_C*k_G*k_X)
FAK=50 - PaxsFAK - PaxFAK
GIT=2 - PaxsGIT - PaxGIT
Raci=-(Rac*Rac_Square - Rac_Square*1 + K0*Rac*alpha_PAK)/Rac_Square
Rhoi=Rac_Square - Rho


R = Rac / Rac_Square
RacRatio = R
RhoRatio = Rho / Rho_Square
PaxRatio = Paxs / Pax_Square
K_is=1.0/((1.0+k_X*PIX+k_G*k_X*k_C*GIT*PIX*Paxtot*PaxRatio)*(1+alpha_R*R)+&
k_G*k_X*GIT*PIX)
K0=alpha_R*K_is*(1+k_X*PIX+k_G*k_X*k_C*Paxtot*GIT*PIX*PaxRatio)
K=R*K0
I_Ks=I_K*(1.0-K_is*(1+alpha_R*R))
P_i=1.0-PaxRatio*(1+k_G*k_X*k_C*GIT*PIX*PAKtot*K_is*(1+alpha_R*R))
Rbar = (Rac + ((K0*Rac*alpha_PAK)/Rac_Square)) / Rac_Square
Q_R = (I_R+I_Ks)*(L_rho**m/(L_rho**m+RhoRatio**m))
Q_rho = I_rho*(L_R**m/(L_R**m +(Rbar)**m))
serine_kinase = (K + PaxsGIT + PaxGIT)
Q_ps = B*(serine_kinase**n)/(L**n+serine_kinase**n)

J_gamma_2_1=(GIT*PAKtot*PIX*k_C*k_G*k_X*((Rac*alpha_R)/Rac_Square + 1))/&
(Pax_Square*(((Rac*alpha_R)/Rac_Square + 1)*(PIX*k_X + (GIT*PIX*Paxs*Paxtot*&
k_C*k_G*k_X)/Pax_Square + 1) + GIT*PIX*k_G*k_X))
J_gamma_1_4=(GIT**2*PIX**2*Paxtot*Pax_Square*Rac*Rac_Square*alpha_R*alpha_PAK*&
k_C*k_G**2*k_X**2)/(Pax_Square*Rac_Square + Pax_Square*Rac*alpha_R + PIX*&
Pax_Square*Rac_Square*k_X + PIX*Pax_Square*Rac*alpha_R*k_X + GIT*PIX*&
Pax_Square*Rac_Square*k_G*k_X + GIT*PIX*Paxs*Paxtot*Rac_Square*k_C*k_G*k_X + &
GIT*PIX*Paxs*Paxtot*Rac*alpha_R*k_C*k_G*k_X)**2
J_gamma_2_4=-(GIT**2*PAKtot*PIX**2*Pax*Paxtot*k_C**2*k_G**2*k_X**2*((Rac*alpha_R)/&
Rac_Square + 1)**2)/(Pax_Square**2*(((Rac*alpha_R)/Rac_Square + 1)*(PIX*k_X + &
(GIT*PIX*Paxs*Paxtot*k_C*k_G*k_X)/Pax_Square + 1) + GIT*PIX*k_G*k_X)**2)
J_gamma_1_6=-(PIX*Pax_Square**2*Rac*Rac_Square*alpha_R*alpha_PAK*k_G*k_X*(PIX*&
k_X + 1))/(Pax_Square*Rac_Square + Pax_Square*Rac*alpha_R + PIX*Pax_Square*&
Rac_Square*k_X + PIX*Pax_Square*Rac*alpha_R*k_X + GIT*PIX*Pax_Square*&
Rac_Square*k_G*k_X + GIT*PIX*Paxs*Paxtot*Rac_Square*k_C*k_G*k_X + GIT*PIX*&
Paxs*Paxtot*Rac*alpha_R*k_C*k_G*k_X)**2
J_gamma_2_6=(PAKtot*PIX*Pax*Pax_Square*k_C*k_G*k_X*(PIX*k_X + 1)*(Rac_Square &
+ Rac*alpha_R)**2)/(Pax_Square*Rac_Square + Pax_Square*Rac*alpha_R + PIX*&
Pax_Square*Rac_Square*k_X + PIX*Pax_Square*Rac*alpha_R*k_X + GIT*PIX*&
Pax_Square*Rac_Square*k_G*k_X + GIT*PIX*Paxs*Paxtot*Rac_Square*k_C*k_G*k_X + &
GIT*PIX*Paxs*Paxtot*Rac*alpha_R*k_C*k_G*k_X)**2
J_gamma_1_10=(alpha_R*alpha_PAK*(PIX*k_X + (GIT*PIX*Paxs*Paxtot*k_C*k_G*k_X)/&
Pax_Square + 1))/(Rac_Square*(((Rac*alpha_R)/Rac_Square + 1)*(PIX*k_X + (GIT*&
PIX*Paxs*Paxtot*k_C*k_G*k_X)/Pax_Square + 1) + GIT*PIX*k_G*k_X)) - (Rac*&
alpha_R**2*alpha_PAK*(PIX*k_X + (GIT*PIX*Paxs*Paxtot*k_C*k_G*k_X)/Pax_Square +&
 1)**2)/(Rac_Square**2*(((Rac*alpha_R)/Rac_Square + 1)*(PIX*k_X + (GIT*PIX*&
Paxs*Paxtot*k_C*k_G*k_X)/Pax_Square + 1) + GIT*PIX*k_G*k_X)**2)
J_gamma_2_10=(GIT**2*PAKtot*PIX**2*Pax*Pax_Square*Rac_Square*alpha_R*k_C*k_G**2*&
k_X**2)/(Pax_Square*Rac_Square + Pax_Square*Rac*alpha_R + PIX*Pax_Square*&
Rac_Square*k_X + PIX*Pax_Square*Rac*alpha_R*k_X + GIT*PIX*Pax_Square*&
Rac_Square*k_G*k_X + GIT*PIX*Paxs*Paxtot*Rac_Square*k_C*k_G*k_X + GIT*PIX*&
Paxs*Paxtot*Rac*alpha_R*k_C*k_G*k_X)**2

f_FAK=PaxFAK*koff1 - FAK*Pax*kon1 - FAK*Paxs*kon1 + PaxsFAK*kfs*koff1
f_PaxFAK=FAK*Pax*kon1 - PaxFAK*koff1
f_Paxs=Pax*Q_ps - Paxs*delta_P + PaxsGIT*koff2 - FAK*Paxs*kon1 - GIT*Paxs*&
kon2 + PaxsFAK*kfs*koff1
f_GIT=PaxGIT*koff2 + PaxsGIT*koff2 - GIT*Pax*kon2 - GIT*Paxs*kon2
f_PaxGIT=GIT*Pax*kon2 - PaxGIT*koff2
f_Raci=Rac*delta_R - Q_R*Raci
f_Rhoi=Rho*delta_rho - Q_rho*Rhoi

F(1)=- PaxFAK*koff1 - (kon1*(PaxFAK + PaxsFAK - cnsrv_2)*(Paxs*Pax_Square + &
PaxFAK*Pax_Square + PaxGIT*Pax_Square + PaxsFAK*Pax_Square + PaxsGIT*&
Pax_Square - Pax_Square*cnsrv_1))/(K_is*PAKtot*PIX*PaxGIT*k_C*k_G*k_X - &
Pax_Square + K_is*PAKtot*PIX*PaxsGIT*k_C*k_G*k_X - K_is*PAKtot*PIX*cnsrv_3*&
k_C*k_G*k_X + K_is*PAKtot*PIX*PaxGIT*R*alpha_R*k_C*k_G*k_X + K_is*PAKtot*PIX*&
PaxsGIT*R*alpha_R*k_C*k_G*k_X - K_is*PAKtot*PIX*R*alpha_R*cnsrv_3*k_C*k_G*k_X)
F(2)=PaxsGIT*koff2 - Paxs*delta_P + (Q_ps*(Paxs*Pax_Square + PaxFAK*&
Pax_Square + PaxGIT*Pax_Square + PaxsFAK*Pax_Square + PaxsGIT*Pax_Square - &
Pax_Square*cnsrv_1))/(K_is*PAKtot*PIX*PaxGIT*k_C*k_G*k_X - Pax_Square + K_is*&
PAKtot*PIX*PaxsGIT*k_C*k_G*k_X - K_is*PAKtot*PIX*cnsrv_3*k_C*k_G*k_X + K_is*&
PAKtot*PIX*PaxGIT*R*alpha_R*k_C*k_G*k_X + K_is*PAKtot*PIX*PaxsGIT*R*alpha_R*&
k_C*k_G*k_X - K_is*PAKtot*PIX*R*alpha_R*cnsrv_3*k_C*k_G*k_X) + PaxsFAK*kfs*&
koff1 + Paxs*kon1*(PaxFAK + PaxsFAK - cnsrv_2) + Paxs*kon2*(PaxGIT + PaxsGIT &
- cnsrv_3)
F(3)=- PaxsFAK*kfs*koff1 - Paxs*kon1*(PaxFAK + PaxsFAK - cnsrv_2)
F(4)=- PaxGIT*koff2 - (kon2*(PaxGIT + PaxsGIT - cnsrv_3)*(Paxs*Pax_Square + &
PaxFAK*Pax_Square + PaxGIT*Pax_Square + PaxsFAK*Pax_Square + PaxsGIT*&
Pax_Square - Pax_Square*cnsrv_1))/(K_is*PAKtot*PIX*PaxGIT*k_C*k_G*k_X - &
Pax_Square + K_is*PAKtot*PIX*PaxsGIT*k_C*k_G*k_X - K_is*PAKtot*PIX*cnsrv_3*&
k_C*k_G*k_X + K_is*PAKtot*PIX*PaxGIT*R*alpha_R*k_C*k_G*k_X + K_is*PAKtot*PIX*&
PaxsGIT*R*alpha_R*k_C*k_G*k_X - K_is*PAKtot*PIX*R*alpha_R*cnsrv_3*k_C*k_G*k_X)
F(5)=- PaxsGIT*koff2 - Paxs*kon2*(PaxGIT + PaxsGIT - cnsrv_3)
F(6)=-(Rac*delta_R + (Q_R*(Rac*Rac_Square - Rac_Square*cnsrv_4 + K0*Rac*&
alpha_PAK))/Rac_Square - (PIX*Pax_Square**2*Rac*Rac_Square*alpha_R*alpha_PAK*&
k_G*k_X*(PIX*k_X + 1)*(PaxGIT*koff2 + PaxsGIT*koff2 + Paxs*kon2*(PaxGIT + &
PaxsGIT - cnsrv_3) + (kon2*(PaxGIT + PaxsGIT - cnsrv_3)*(Paxs*Pax_Square + &
PaxFAK*Pax_Square + PaxGIT*Pax_Square + PaxsFAK*Pax_Square + PaxsGIT*&
Pax_Square - Pax_Square*cnsrv_1))/(K_is*PAKtot*PIX*PaxGIT*k_C*k_G*k_X - &
Pax_Square + K_is*PAKtot*PIX*PaxsGIT*k_C*k_G*k_X - K_is*PAKtot*PIX*cnsrv_3*&
k_C*k_G*k_X + K_is*PAKtot*PIX*PaxGIT*R*alpha_R*k_C*k_G*k_X + K_is*PAKtot*PIX*&
PaxsGIT*R*alpha_R*k_C*k_G*k_X - K_is*PAKtot*PIX*R*alpha_R*cnsrv_3*k_C*k_G*&
k_X)))/(Pax_Square*Rac_Square + Pax_Square*Rac*alpha_R + PIX*Pax_Square*&
Rac_Square*k_X + PIX*Pax_Square*Rac*alpha_R*k_X - PIX*Pax_Square*Rac_Square*&
k_G*k_X*(PaxGIT + PaxsGIT - cnsrv_3) - PIX*Paxs*Paxtot*Rac_Square*k_C*k_G*&
k_X*(PaxGIT + PaxsGIT - cnsrv_3) - PIX*Paxs*Paxtot*Rac*alpha_R*k_C*k_G*k_X*&
(PaxGIT + PaxsGIT - cnsrv_3))**2 + (PIX**2*Paxtot*Pax_Square*Rac*Rac_Square*&
alpha_R*alpha_PAK*k_C*k_G**2*k_X**2*(PaxGIT + PaxsGIT - cnsrv_3)**2*(PaxsGIT*&
koff2 - Paxs*delta_P + (Q_ps*(Paxs*Pax_Square + PaxFAK*Pax_Square + PaxGIT*&
Pax_Square + PaxsFAK*Pax_Square + PaxsGIT*Pax_Square - Pax_Square*cnsrv_1))/&
(K_is*PAKtot*PIX*PaxGIT*k_C*k_G*k_X - Pax_Square + K_is*PAKtot*PIX*PaxsGIT*&
k_C*k_G*k_X - K_is*PAKtot*PIX*cnsrv_3*k_C*k_G*k_X + K_is*PAKtot*PIX*PaxGIT*R*&
alpha_R*k_C*k_G*k_X + K_is*PAKtot*PIX*PaxsGIT*R*alpha_R*k_C*k_G*k_X - K_is*&
PAKtot*PIX*R*alpha_R*cnsrv_3*k_C*k_G*k_X) + PaxsFAK*kfs*koff1 + Paxs*kon1*&
(PaxFAK + PaxsFAK - cnsrv_2) + Paxs*kon2*(PaxGIT + PaxsGIT - cnsrv_3)))/&
(Pax_Square*Rac_Square + Pax_Square*Rac*alpha_R + PIX*Pax_Square*Rac_Square*&
k_X + PIX*Pax_Square*Rac*alpha_R*k_X - PIX*Pax_Square*Rac_Square*k_G*k_X*&
(PaxGIT + PaxsGIT - cnsrv_3) - PIX*Paxs*Paxtot*Rac_Square*k_C*k_G*k_X*&
(PaxGIT + PaxsGIT - cnsrv_3) - PIX*Paxs*Paxtot*Rac*alpha_R*k_C*k_G*k_X*&
(PaxGIT + PaxsGIT - cnsrv_3))**2)/((alpha_R*alpha_PAK*(PIX*k_X - (PIX*Paxs*&
Paxtot*k_C*k_G*k_X*(PaxGIT + PaxsGIT - cnsrv_3))/Pax_Square + 1))/&
(Rac_Square*(((Rac*alpha_R)/Rac_Square + 1)*(PIX*k_X - (PIX*Paxs*Paxtot*k_C*&
k_G*k_X*(PaxGIT + PaxsGIT - cnsrv_3))/Pax_Square + 1) - PIX*k_G*k_X*(PaxGIT +&
 PaxsGIT - cnsrv_3))) - (Rac*alpha_R**2*alpha_PAK*(PIX*k_X - (PIX*Paxs*Paxtot*&
k_C*k_G*k_X*(PaxGIT + PaxsGIT - cnsrv_3))/Pax_Square + 1)**2)/(Rac_Square**2*&
(((Rac*alpha_R)/Rac_Square + 1)*(PIX*k_X - (PIX*Paxs*Paxtot*k_C*k_G*k_X*&
(PaxGIT + PaxsGIT - cnsrv_3))/Pax_Square + 1) - PIX*k_G*k_X*(PaxGIT + &
PaxsGIT - cnsrv_3))**2) + 1)
F(7)=- Rho*delta_rho - Q_rho*(Rho - cnsrv_5)
END SUBROUTINE FUNC

SUBROUTINE STPNT(NDIM,U,PAR,T)
IMPLICIT NONE
INTEGER NDIM
DOUBLE PRECISION U(NDIM), PAR(*), T

PAR(1)=1.000000000 !kon1
PAR(2)=1.000000000 !koff1
PAR(3)=1.000000000 !kon2
PAR(4)=1.000000000 !koff2
PAR(5)=50.000000000 !kfs
PAR(6)=0.750000000 !L
PAR(7)=4.000000000 !n
PAR(8)=2.000000000 !B
PAR(9)=0.016000000 !I_rho
PAR(10)=0.340000000 !L_rho
PAR(12)=0.016000000 !delta_rho
PAR(13)=0.340000000 !L_R
PAR(14)=0.003000000 !I_R
PAR(15)=0.025000000 !delta_R
PAR(16)=15.000000000 !alpha_R
PAR(17)=0.000400000 !delta_P
PAR(18)=0.009000000 !I_K
PAR(19)=41.700000000 !k_X
PAR(20)=5.710000000 !k_G
PAR(21)=5.000000000 !k_C
PAR(22)=0.069000000 !PIX
PAR(23)=2.300000000 !Paxtot
PAR(24)=4.000000000 !m
PAR(25)=0.300000000 !alpha_PAK
PAR(26)=1.000000000 !Rho_Square
PAR(27)=1.000000000 !Rac_Square
PAR(28)=1.000000000 !Pax_Square
PAR(29)=2.250000000 !PAKtot

U(1)= 0.00325719587762
U(2)=   0.280146333675
U(3)=   0.278567291736
U(4)=0.000102347407723
U(5)=   0.437656211919
U(6)=   0.235501772759
U(7)=   0.192084565275
END SUBROUTINE STPNT

SUBROUTINE BCND
END SUBROUTINE BCND

SUBROUTINE ICND
END SUBROUTINE ICND

SUBROUTINE FOPT
END SUBROUTINE FOPT

SUBROUTINE PVLS
END SUBROUTINE PVLS

