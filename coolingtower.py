#Cooling Tower Model
#Shreyas Sudhakar
#2018

import CoolProp.CoolProp as CP
import numpy as np
from scipy.optimize import root

#Define known variables
load = 1000 #Watts, cooling load from business building
height = 3 #Meters, cooling tower height
A_c = 2 #Square meters, cooling tower cross sectional area
P_0 = 101325 #Pa, atmospheric pressure
Dab = 0.26e-4 #m^2/s Water air diffusivity atmospheric

dx = 2 #Distance step

#Boundary Conditions
T_a_0 = 295.15 #Kelvin, atmopheric air temp AKA T_a(X=0)
W_a_0 = 0.005 #Humidity Ratio atmospheric AKA W_a(X=0)
t_F_n = 312 #Kelvin, hot water temperature AKA T_F(X=n)
m_dot_W_n = 1 #L/min, hot water mass flow rate AKA m_dot_W(X=n)

x_0 = 0 #Meters, base of cooling tower (X=0)
x_n = 0 + height #Meters, top of cooling tower (X=n)

#estimate water quantity
v_total = A_c*height #m^3, Total volume of the cooling tower
P_v = CP.HAPropsSI('P_w','W',W_a_0,'T',T_a_0,'P',P_0)
v_w = (P_v/P_0)*v_total #m^3, water volume of cooling tower computed using partial pressures


#Droplet properties
d_p_0 = 5E-6 #m, droplet starting diameter. Assumption based on https://www.ncbi.nlm.nih.gov/pubmed/1240052
v_p_0 = (4/3)*np.pi*(d_p_0/2)**3 #m^3, initial droplet volume
A_p_0 = 4*np.pi*(d_p_0/2)**2 #m^2, initial droplet surface area
N_p = v_w/v_p_0 #Number of particles, assuming that this is a constant
A_s_0 = A_p_0*N_p #Surface area of all the droplates

#A)
#Calculate the inlet air flow rate
#Rules of thumb states air flow rates between 6344 - 8784 kg/h.
m_dot_A_0 = np.mean([6344,8784]) #kg/h, Using the average expected air flow rate
rhoA = CP.PropsSI("D", "T", T_a_0, "P", P_0, "Air") #kg/m^3, Compute air density at atmospheric conditions
rhoinf = CP.PropsSI("D", "T", T_a_0, "Q", 1, "Water")#kg/m^3, Compute water density at atmospheric conditions
Q_A_0 = m_dot_A_0/rhoA #m^3/h, air volumetric flow rate
Q_A_0 = Q_A_0/3600 #m^3/s, air volumetric flow rate
velocity = Q_A_0/A_c #m/s, air velocity


viscosityAir = CP.PropsSI("V", "T", T_a_0, "P", P_0, "Air") #Pa-s, viscosity of air
v_a = v_total-v_w #m^3, air volume
m_a = v_a*rhoA #kg, air mass

c_p = CP.PropsSI("C", "T", T_a_0, "P", P_0, "Air") #J/kg/K Cp mass specific constant pressure specific heat, assuming this is constant

#outlet water temperature
#Per "Cooling Tower Rules of Thumb, Countercurrent-induced draft towers are the most common in process industries. They are able to cool water within 2°F of the wet bulb.
t_outlet_water = CP.HAPropsSI('WetBulb','T',T_a_0,'P',101325,'W',W_a_0) + 2 #Kelvin, Outlet water temperature
print(t_outlet_water)


def properties(inputs):
    m_dot = inputs[0:int(x_n/dx)-1]
    T = inputs[int(x_n/dx):]

    errors = np.ones((4, int(x_n/dx)))
    for distance in range(0, len(m_dot)-1):
        reynolds = velocity * dx/viscosityAir
        Sc = viscosityAir/(rhoA*Dab)
        Sh = 0.43*reynolds**0.58*Sc**0.4
        hm = Dab*Sh/l
        rhoFs = CP.PropsSI("D", "T", T[int(distance)], "Q", 1, "Water")

        #State Calculations
        errors[0][int(distance)] = (m_dot[int(distance)+1] - m_dot[int(distance)])/dx - hm*A_p_0*(rhoFs - rhoFinf)*A_c/0.622/velocity
        errors[1][int(distance)] = (m_dot[int(distance)+1] - m_dot[int(distance)])/dx - m_a*W_a_0 #unsure if this is implemented correctly
        errors[2][int(distance)] = ((m_dot[int(distance)+1] - m_dot[int(distance)])/dx)*c_p*(T_a[int(distance)+1] - T_a[int(distance)]) - hm*A_p_0*(T_a_0-T_f[int(distance)])*A_c/0.622/v_p_0
        errors[3][int(distance)] = m_a*c_p*(T_a[int(distance)+1] - T_a[int(distance)])-((m_dot[int(distance)+1] - m_dot[int(distance)])/dx)*c_p*(T_a[int(distance)+1]

        #Boundary Conditions
        errors[0][-1] = m_dot[-1] - m_dot_W_n
        errors[1][-1] = m_a[-1] - (v_total - (W_a_0*(v_total))
        errors[2][-1] = T_f[-1] - t_outlet_water
        errors[3][-1] = T_a[-1] - T_a_0
    return np.reshape(errors, 4*int(x_n/dx)) #need to be reshaped as errors must be returned as a long array
