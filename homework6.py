import compressibleflows.isentropic as isen
import compressibleflows.shockjump as sj
import compressibleflows.prandltmeyer as prandlt
from compressibleflows import obshock
import numpy as np 

# Problem 1 Anderson 9.12
M1 = 3
T1 = 285 # K
P1 = 1 # atm
theta = 30.6
P01 = isen.P0_P(M1)*P1 
T01 = isen.T0_T(M1)*T1 
beta,M1n,M2,theta_max = obshock(M1,theta)
print(f"Beta={beta:0.2f}, M1n={M1n:0.2f}, M2 ={M2:0.2f},theta_max={theta_max:0.2f}") # After oblique shock
P2_P1 = sj.P2_P1(M1n)
T2_T1 = sj.T2_T1(M1n)
print(f"P2_P1={P2_P1:0.2f} T2_T1={T2_T1:0.2f}")
P2 = P2_P1*P1
T2 = T2_T1*T1
P02 = P01*sj.P02_P01(M1n)
T02 = T01
print(f"P2={P2:0.2f} T2={T2:0.2f}")
# Prandlt Meyer
nu_m2,_ = prandlt.prandlt(M2)
nu_m3 = theta + nu_m2
M3 = prandlt.prandltM(nu_m3)
print(f"M3={M3:0.2f}")
P0_P = isen.P0_P(M3)
T0_T = isen.T0_T(M3)
T3 = T02/T0_T
P3 = P02/P0_P
print(f"P3={P3:0.2f} T3={T3:0.2f}")

# Problem 2
M1 = 2
M3 = 4 
nu_m2 = (prandlt.prandlt(M1)[0]+prandlt.prandlt(M3)[0])/2 
theta = nu_m2 - prandlt.prandlt(M1)[0]
print(f"Theta={theta}")

# Problem 3 
M1 = 2.6
alpha = 30
## Suction side
P01_P1 = isen.P0_P(M1)
nu_m2 = alpha + prandlt.prandlt(M1)[0]
M2 = prandlt.prandltM(nu_m2)
print(f"nu_m2={nu_m2:0.2f} M2={M2:0.2f}")
P02_P2 = isen.P0_P(M2)
P2_P1 = P01_P1/P02_P2
print(f"P2P1={P2_P1:0.2f}")
## Pressure side
beta,M1n,M3,theta_max = obshock(M1,alpha)
print(f"Beta={beta:0.2f}, M1n={M1n:0.2f}, M3 ={M3:0.2f},theta_max={theta_max:0.2f}") # After oblique shock
P3_P1 = sj.P2_P1(M1n)
print(f"P3P1={P3_P1:0.2f}")
L_prime = (P3_P1-P2_P1)/(0.5*1.4)
Cd = L_prime*np.sin(np.radians(alpha))
Cl = L_prime*np.cos(np.radians(alpha))
print(f"Cl={Cl}  Cd={Cd}")