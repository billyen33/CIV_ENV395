#Part 1
import math
import pandas as pd

#generate time array
t = []
for ii in range(50):
    t.append(ii)
#print(t)
#initialize constants
R = 0.082057 #(L*atm)/(mol*K)
Q = 25 #m^3/hr
V = 350
depth = 12
SA = V/depth #m^2
Cin = 50 #mg/m^3

U_denver = 5 * 3600
U_pho = 2* 3600
U_chi = 8* 3600

Temp_d = ((45-32)*5)/9 + 273.15
Temp_p = ((78-32)*5)/9 + 273.15
Temp_c = ((57-32)*5)/9 + 273.15

kd_d = 0.025
kd_p = 0.031
kd_c = 0.034

Da_d = 10.81
Da_p = 4.64
Da_c = 3.5

Dw_d = 5.27*10**(-6)
Dw_p = 8.04*10**(-6)
Dw_c = 1*10**(-5)

Kh_d = 23.5
Kh_p = 12.4
Kh_c = 18.9

Kl_d = ((Dw_d*0.0014*(U_denver)**2)/(2.6*10**(-5)))**(1/3)
Kl_p = ((Dw_p*0.0014*(U_pho)**2)/(2.6*10**(-5)))**(1/3)
Kl_c = ((Dw_c*0.0014*(U_chi)**2)/(2.6*10**(-5)))**(1/3)

Kg_d = ((Da_d/0.26)*7*U_denver)**(1/2)
Kg_p = ((Da_p/0.26)*7*U_pho)**(1/2)
Kg_c = ((Da_c/0.26)*7*U_chi)**(1/2)

Kgl_d = 1/((1/Kl_d)+(Kh_d*R*Temp_d)/Kg_d)
Kgl_p = 1/((1/Kl_p)+(Kh_d*R*Temp_p)/Kg_p)
Kgl_c = 1/((1/Kl_c)+(Kh_c*R*Temp_c)/Kg_c)

#Equation time!!!
S = (Q*Cin)/V #defined source term
Lc = ((Q+Kgl_c*SA)/V)-kd_c
Ld = ((Q+Kgl_d*SA)/V)-kd_d
Lp = ((Q+Kgl_p*SA)/V)-kd_p
#Chicago
C_chicago = []
#generate C with respect to time
for item in t:
    C_chicago.append(S/Lc - (S/Lc - Cin)*math.exp(-Lc*item))
#Denver
C_denver = []
#generate C with respect to time
for item in t:
    C_denver.append(S/Lc - (S/Ld - Cin)*math.exp(-Ld*item))
#Pheonix
C_phoenix = []
#generate C with respect to time
for item in t:
    C_phoenix.append(S/Lc - (S/Lp - Cin)*math.exp(-Lp*item))

#graph
from matplotlib import pyplot as plt
fig, ax = plt.subplots()
ax.plot(t, C_chicago, label = 'Chicago')
ax.plot(t, C_denver, label = "Denver")
ax.plot(t, C_phoenix, label = 'Phoenix')
ax.set_xlabel('Time (hours)')
ax.set_ylabel('Concentration of MTBE (mg/m^3)')
ax.set_title('Concentration of MTBE (mg/m^3) vs Time (hours) in Three Cities')
ax.legend()
plt.show()