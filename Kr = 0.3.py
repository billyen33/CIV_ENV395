from Case_Study_2 import Y
from matplotlib import pyplot as plt
import math
import numpy as np
Kr = 0.3/24
U_ET = 1999.25
Kbod_ev = 9.05/24 #0.4605/24
Lo_ET = 22347.949023498655
D_not_ET = 5884.675610433376
DOir = 9000
Xc_ET = (U_ET/(Kr-Kbod_ev))*np.log(Kr/Kbod_ev)*(1-(D_not_ET)*((Kr-Kbod_ev)/(Kbod_ev*Lo_ET)))

distance_ET_bad = []
DOx_ET = []
for item in range(int(Xc_ET)+10001):
    distance_ET_bad.append(item)
    Dx_ET = (((Kbod_ev*Lo_ET)/(Kr-Kbod_ev))*(math.exp((-Kbod_ev*item)/U_ET)-math.exp((-Kr*item)/U_ET))+D_not_ET*math.exp((-Kr*item)/U_ET))
    DOx_ET.append(DOir - Dx_ET)
D = (Kbod_ev*Lo_ET*math.exp(-Kbod_ev*Xc_ET/U_ET))/(Kr)
print('theoretical DOmin ' + str(DOir - D))
print('min of graph ' + str(min(DOx_ET)))
print('Index: '+ str(DOx_ET.index(min(DOx_ET))))
print(Xc_ET)
fig3, ax3 = plt.subplots()
#ax3.plot(distance_ET_bad, DO_limit, linestyle = 'dashed', color='red', label = 'Safety Level')
ax3.plot(distance_ET_bad, DOx_ET, label = 'Sag Curve')
ax3.axvline(x = Xc_ET, ymin = 0, ymax = 9, linestyle = '--', color = 'gray', alpha = 0.8, label='Xc of EvansTown')
ax3.set_xlabel('Distance (km)')
ax3.set_ylabel('Concentration of DO (mg/L)')
ax3.set_title('DO Sag Curve ET')
ax3.legend()

plt.show()