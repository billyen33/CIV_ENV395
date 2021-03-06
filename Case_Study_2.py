#Part 1
import math
import numpy as np
#initialize constants
R = 0.082057 #(L*atm)/(mol*K)
Q = 25 #m^3/hr
V = 350 #m^3
depth = 12 #m
SA = V/depth #m^2
Cin = 50 #mg/m^3
theta = V/Q #hour

#generate time list
t1 = []
for ii in range(160): #part 1
    t1.append(ii/10)
t2 = []
for ii in range(100): #part 2
    t2.append(ii)

U_denver = 5 * 3600 #in m/hours
U_pho = 2 * 3600
U_chi = 8 * 3600

Temp_d = ((45-32)*5)/9 + 273.15 #convert F to K
Temp_p = ((78-32)*5)/9 + 273.15
Temp_c = ((57-32)*5)/9 + 273.15

kd_d = 0.025
kd_p = 0.031
kd_c = 0.034

Da_d = 10.81
Da_p = 4.64
Da_c = 3.5

Dw_d = 5.27*(10**(-6))
Dw_p = 8.04*(10**(-6))
Dw_c = 1*(10**(-5))

Kh_d = 23.5
Kh_p = 12.4
Kh_c = 18.9

Kl_d = ((Dw_d*0.0014*((U_denver)**2))/(2.6*(10**(-5))))**(1/3)
Kl_p = ((Dw_p*0.0014*((U_pho)**2))/(2.6*(10**(-5))))**(1/3)
Kl_c = ((Dw_c*0.0014*((U_chi)**2))/(2.6*(10**(-5))))**(1/3)

Kg_d = ((Da_d/0.26)*7*U_denver)**(1/2)
Kg_p = ((Da_p/0.26)*7*U_pho)**(1/2)
Kg_c = ((Da_c/0.26)*7*U_chi)**(1/2)
Kgl_d = 1/((1/Kl_d)+((Kh_d*R*Temp_d)/Kg_d))
Kgl_p = 1/((1/Kl_p)+((Kh_d*R*Temp_p)/Kg_p))
Kgl_c = 1/((1/Kl_c)+((Kh_c*R*Temp_c)/Kg_c))
'''print(Kgl_d)
print(Kgl_p)
print(Kgl_c)'''
#Equation time!!!
S = (Q*Cin)/V #defined source term
Lc = ((Q+Kgl_c*SA)/V)+kd_c #loss terms defined here
Ld = ((Q+Kgl_d*SA)/V)+kd_d
Lp = ((Q+Kgl_p*SA)/V)+kd_p
#Chicago
C_chicago = []
#generate C with respect to time
for item in t1:
    C_chicago.append(S/Lc - (S/Lc - Cin)*math.exp(-Lc*item))

#Denver
C_denver = []
#generate C with respect to time
for item in t1:
    C_denver.append(S/Ld - (S/Ld - Cin)*math.exp(-Ld*item))

#Pheonix
C_phoenix = []
#generate C with respect to time
for item in t1:
    C_phoenix.append(S/Lp - (S/Lp - Cin)*math.exp(-Lp*item))
#print('Residence time: ' + str(theta) + 'hrs')
#print('Phoenix: '+ str(S/Lp - (S/Lp - Cin)*math.exp(-Lp*(theta))))
#print('Chicago: '+ str(S/Lc - (S/Lc - Cin)*math.exp(-Lc*(theta))))
#print('Not our Chicago: ' + str(S/Lc))
#print('Denver: '+ str(S/Ld - (S/Ld - Cin)*math.exp(-Ld*(theta))))
#graph
safe_MTBE = []
for item in t1:
    safe_MTBE.append(10)
from matplotlib import pyplot as plt
fig, ax = plt.subplots()
ax.plot(t1, safe_MTBE, linestyle = 'dashed', color='red', label = 'Safety Level')
ax.plot(t1, C_chicago, label = 'Chicago')
ax.plot(t1, C_denver, label = "Denver")
ax.plot(t1, C_phoenix, label = 'Phoenix')
ax.axvline(x = V/Q, ymin = 0, ymax = 9, linestyle = ':', color = 'gray', alpha = 0.8, label='Residence Time')
ax.set_xlabel('Time (hours)')
ax.set_ylabel('Concentration of MTBE (mg/m^3)')
ax.set_title('Concentration of MTBE (mg/m^3) vs Time (hours) in Three Cities')
ax.legend()

#Part 2
DOmin = 4 * 1000 #g/m^3 * 1000 mg/g
DOw = 1 * 1000 #g/m^3 * 1000 mg/g
DOs = 9 * 1000 #g/m^3 * 1000 mg/g
Y = 2.7273
Kglo_d = 0.25 #m/hr
Kglo_p = 1.2
Kglo_c = 1
DOt_d = []
DOt_p = []
DOt_c = []
# Solve for DO(t)
for item in t2:
    #Cmtbe at each city
    Cmtbe_d = S/Ld - (S/Ld - Cin)*math.exp(-Ld*(item))
    Cmtbe_p = S/Lp - (S/Lp - Cin)*math.exp(-Lp*(item))
    Cmtbe_c = S/Lc - (S/Lc - Cin)*math.exp(-Lc*(item))
    #S terms
    '''Soxygen_d = (Q*DOw/V) - (X*kd_d*Cmtbe_d)
    Soxygen_p = (Q*DOw/V) - (X*kd_p*Cmtbe_p)
    Soxygen_c = (Q*DOw/V) - (X*kd_c*Cmtbe_c)'''

    '''Soxygen_d = (Q*DOw/V) - (X*kd_d*Cmtbe_d) + (Kglo_d*DOs*SA)/V
    Soxygen_p = (Q*DOw/V) - (X*kd_p*Cmtbe_p) + (Kglo_p*DOs*SA)/V
    Soxygen_c = (Q*DOw/V) - (X*kd_c*Cmtbe_c) + (Kglo_c*DOs*SA)/V'''
    Soxygen_d = (Q*DOw/V) - (Y*kd_d*Cmtbe_d) + (Kglo_d*DOs*SA)/V
    Soxygen_p = (Q*DOw/V) - (Y*kd_p*Cmtbe_p) + (Kglo_p*DOs*SA)/V
    Soxygen_c = (Q*DOw/V) - (Y*kd_c*Cmtbe_c) + (Kglo_c*DOs*SA)/V
    #L terms
    '''Loxygen_d = (-Kglo_d*SA+Q)/V
    Loxygen_p = (-Kglo_p*SA+Q)/V
    Loxygen_c = (-Kglo_c*SA+Q)/V'''

    '''Loxygen_d = Q/V
    Loxygen_p = Q/V
    Loxygen_c = Q/V'''
    Loxygen_d = (Kglo_d*SA+Q)/V
    Loxygen_p = (Kglo_p*SA+Q)/V
    Loxygen_c = (Kglo_c*SA+Q)/V
    #DOu terms
    DOu_d = Soxygen_d/Loxygen_d
    DOu_p = Soxygen_p/Loxygen_p
    DOu_c = Soxygen_c/Loxygen_c
    #solve for DO(t)
    DOt_d.append((DOu_d - (DOu_d - DOs)*math.exp(-Loxygen_d*item))/1000) #add in terms of g/m^3
    DOt_p.append((DOu_p - (DOu_p - DOs)*math.exp(-Loxygen_p*item))/1000)
    DOt_c.append((DOu_c - (DOu_c - DOs)*math.exp(-Loxygen_c*item))/1000)
t = 24*5 #set time to large number to get C infinity so we can calc BODu of pond
Cmtbe_c = S/Lc - (S/Lc - Cin)*math.exp(-Lc*(t))
Soxygen_c = (Q*DOw/V) - (Y*kd_c*Cmtbe_c) + (Kglo_c*DOs*SA)/V
Loxygen_c = (Kglo_c*SA+Q)/V
DOu_c = Soxygen_c/Loxygen_c

#print(DOt_c[50])
#print(DOt_c)
#graph
line = []
for item in t2:
    line.append(4)
fig1, ax1 = plt.subplots()
ax1.plot(t2, line, linestyle = 'dashed', color='red', label = 'Safety Level')
ax1.plot(t2, DOt_c, label = 'Chicago')
ax1.plot(t2, DOt_d, label = "Denver")
ax1.plot(t2, DOt_p, label = 'Phoenix')
ax1.set_xlabel('Time (hours)')
ax1.set_ylabel('Concentration of DO (mg/L)')
ax1.set_title('Concentration of DO (mg/L) vs Time (hours) in Three Cities')
ax1.legend()


#Part 3
Kbod_c = 0.2155/24 #1/day * 1 day/ 24hr
Kbod_tv = 0.2197/24 #1/day
Kbod_ev = 0.4605/24
Kr = 0.3/24 #1/day
DOir = 9 * 1000 #mg/L * 1000L/1m^3
DOiw = 9 *1000
DOitv = 9 * 1000
DOiet = 9 * 1000

DOur = 7 * 1000#mg/L* 1000L/1m^3
DOuw = 4.3 * 1000
DOutv = 6 * 1000
DOuet = 5 * 1000

Fr = Ffg = 1
Fw = Ftv = Fet = 0.05

Qr = 9*3600 #m^3/s * 3600 s/hr
Qw = 1.2*3600
Qtv = 0.9*3600
Qet = 1.1*3600
Av = 20 #m^2
U_c = (Qr+Qw+Q)/Av
D_not = DOir - (DOir*Qr + DOiw*Qw + DOt_c[-1]*Q)/(Q+Qr+Qw)
Cmtbe_theta = S/Lc - (S/Lc - Cin)*math.exp(-Lc*(theta))
#From Chicago to TechVille
Lo = (((DOir-DOur)/Fr)*Qr+((DOiw-DOuw)/Fw)*Qw+Y*Cmtbe_theta*Q)/(Qr+Qw+Q)
Xc_c = (U_c/(Kr-Kbod_c))*np.log(Kr/Kbod_c)*(1-(D_not*(Kr-Kbod_c))/(Kbod_c*Lo))
distance_c = []
DOx_c = []
for item in range(int(Xc_c)+6001):
    distance_c.append(item)
    Dx_c = (((Kbod_c*Lo)/(Kr-Kbod_c))*(math.exp((-Kbod_c*item)/U_c)-math.exp((-Kr*item)/U_c))+D_not*math.exp((-Kr*item)/U_c))
    DOx_c.append(DOir - Dx_c)
'''print('Crit Distance Chicago')
print(DOx_c.index(min(DOx_c)))
print(Xc_c)'''

#From TechVille to EvansTown
U_t = (Qr+Qw+Q)/Av
DOi_into_tv = DOx_c[int(Xc_c)+6000]
D_not_tv = DOir - (DOi_into_tv*(Qw+Q+Qr)+DOitv*Qtv)/(Qw+Q+Qr+Qtv)
F_total_chicago = (Fr*Qr+Ffg*Q+Fw*Qw)/(Qr+Q+Qw)
DOu_fromChicago = DOi_into_tv - Lo*F_total_chicago
Lo_tv = (((DOi_into_tv-DOu_fromChicago)/F_total_chicago)*(Q+Qw+Qr)+((DOitv-DOutv)/Ftv)*Qtv)/(Qr+Qw+Q+Qtv)
Xc_tv = (U_t/(Kr-Kbod_tv))*np.log(Kr/Kbod_tv)*(1-(D_not_tv*(Kr-Kbod_tv))/(Kbod_tv*Lo_tv))
distance_tv_bad = []
DOx_tv = []

for item in range(int(Xc_tv)-7999):
    distance_tv_bad.append(item)
    Dx_tv = (((Kbod_tv*Lo_tv)/(Kr-Kbod_tv))*(math.exp((-Kbod_tv*item)/U_t)-math.exp((-Kr*item)/U_t))+D_not_tv*math.exp((-Kr*item)/U_t))
    DOx_tv.append(DOir - Dx_tv)
'''print('Crit Distance TV')
print(DOx_tv.index(min(DOx_tv)))
print(Xc_tv)'''
distance_tv = []
for item in distance_tv_bad:
    distance_tv.append(item+Xc_c+6000)

#From EvansTown onward
U_ET = (Qr+Qw+Q+Qtv)/Av
DOi_into_ET = DOx_tv[int(Xc_tv)-8000]
print(int(Xc_c) + 6000 + int(Xc_tv)-8000)
print(len(DOx_tv))
print(DOi_into_ET)
D_not_ET = DOir - (DOi_into_ET*(Qw+Q+Qr+Qtv)+DOiet*Qet)/(Qw+Q+Qr+Qtv+Qet)
F_total_TV = (Fr*Qr+Ffg*Q+Fw*Qw+Ftv*Qtv)/(Qr+Q+Qw+Qtv)
DOu_fromTechVille = DOi_into_ET - Lo_tv*F_total_TV
Lo_ET = (((DOi_into_ET-DOu_fromTechVille)/F_total_TV)*(Q+Qw+Qr+Qtv)+((DOiet-DOuet)/Fet)*Qet)/(Qr+Qw+Q+Qtv+Qet)
#Lo_ET = 22347
Xc_ET = (U_ET/(Kr-Kbod_ev))*np.log(Kr/Kbod_ev)*(1-(D_not_ET*(Kr-Kbod_ev))/(Kbod_ev*Lo_ET))

distance_ET_bad = []
DOx_ET = []
for item in range(int(Xc_ET)+10001):
    distance_ET_bad.append(item)
    Dx_ET = (((Kbod_ev*Lo_ET)/(Kr-Kbod_ev))*(math.exp((-Kbod_ev*item)/U_ET)-math.exp((-Kr*item)/U_ET))+D_not_ET*math.exp((-Kr*item)/U_ET))
    DOx_ET.append(DOir - Dx_ET)
distance_ET = []
for item in distance_ET_bad:
    distance_ET.append(item+distance_tv[-1])
'''print('Crit Distance ET')
print(DOx_ET.index(min(DOx_ET)))
print(Xc_ET)'''
#total distance list generated here
total_d_bad = distance_c + distance_tv + distance_ET
total_d = []
for item in total_d_bad:
    total_d.append(item/1000)
total_DOx = DOx_c + DOx_tv +DOx_ET
#unit conversion
for ii in range(len(total_DOx)):
    total_DOx[ii] = total_DOx[ii]/1000
DO_limit = []
for item in total_d:
    DO_limit.append(4)

'''print('Do of Chicago ' + str(D_not))
print('Lo of ET ' + str(Lo_ET))
print('Do of EvansTown ' + str(D_not_ET))'''

fig2, ax2 = plt.subplots()
ax2.plot(total_d, DO_limit, linestyle = 'dashed', color='red', label = 'Safety Level')
ax2.plot(total_d, total_DOx, label = 'Sag Curve')
ax2.axvline(x = Xc_c/1000, ymin = 0, ymax = 9, linestyle = ':', color = 'gray', alpha = 0.8, label='Xc of Chicago')
ax2.axvline(x = Xc_c/1000 + 6 + Xc_tv/1000, ymin = 0, ymax = 9, linestyle = '-.', color = 'gray', alpha = 0.8, label='Xc of TechVille')
ax2.axvline(x = Xc_c/1000 + 6 + Xc_tv/1000 - 8 + Xc_ET/1000, ymin = 0, ymax = 9, linestyle = '--', color = 'gray', alpha = 0.8, label='Xc of EvansTown')
#ax2.axvline(x = DOx_c.index(min(DOx_c))/1000 + 6 + DOx_tv.index(min(DOx_tv))/1000 - 8 + DOx_ET.index(min(DOx_ET))/1000, ymin = 0, ymax = 9, linestyle = '--', color = 'gray', alpha = 0.8, label='Xc of EvansTown')
ax2.set_xlabel('Distance (km)')
ax2.set_ylabel('Concentration of DO (mg/L)')
ax2.set_title('DO Sag Curve')
ax2.legend()
'''
fig3, ax3 = plt.subplots()
ax3.plot(total_d, DO_limit, linestyle = 'dashed', color='red', label = 'Safety Level')
ax3.plot(distance_tv_bad, DOx_tv, label = 'Sag Curve')
ax3.axvline(x = Xc_tv, ymin = 0, ymax = 9, linestyle = '--', color = 'gray', alpha = 0.8, label='Xc of EvansTown')
ax3.set_xlabel('Distance (km)')
ax3.set_ylabel('Concentration of DO (mg/L)')
ax3.set_title('DO Sag Curve ET')
ax3.legend()'''

plt.show()