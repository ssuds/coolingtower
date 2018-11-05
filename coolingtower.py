#Cooling Tower Model
#Shreyas Sudhakar
#2018

import CoolProp.CoolProp as CP
import numpy as np
from scipy.optimize import root

def Hf(T):
    if  (T < 273.15 or T > 645 ): #solution moved out of coolprop bounds
        return -1e6
    else:
        return CP.PropsSI('H', 'T', T, 'P', 101325, 'Water') #Computing enthalpy of coolant fluid

def Hfg(T):
    if  (T < 273.15 or T > 645 ): #solution moved out of coolprop bounds
        return -1e6
    else:
        return CP.PropsSI('H', 'T', T, 'Q', 1, 'Water') - CP.PropsSI('H', 'T', T, 'Q', 0, 'Water') #Computing enthalpy of vaporization

def Psat(T):
    if  (T < 273.15): #solution moved out of coolprop bounds
        return -1e6
    else:
        return CP.PropsSI('P', 'T', T, 'Q', 1, 'Water') #computing saturation pressure

def rhoG(T):
    if  (T < 273.15 or T > 645 ): #solution moved out of coolprop bounds
        return -1e6
    else:
        return CP.PropsSI('D', 'T', T, 'Q', 1, 'Water') #compute density

def wSat(T):
    return CP.HAPropsSI('W', 'T', T, 'P', 101325, 'R', 1) #compute humidity ratio for saturated air

def Rh(T, w):
    if  (T < 273.15 or T > 645 ): #solution moved out of coolprop bounds
        return -1e6
    else:
        return CP.HAPropsSI('R', 'T', T, 'P', 101325, 'W', w) #compute relative humidity given temperature and humidity ratio

def Ha(T, w):
    if  (T < 273.15 or T > 645 ): #solution moved out of coolprop bounds
        return -1e6
    else:
        return CP.HAPropsSI('H', 'T', T, 'P', 101325, 'W', w) #compute enthalpy of air given temperature and humidity ratio




#Define known variables
load = 1000 #Watts, cooling load from business building
height = 3 #Meters, cooling tower height (aka towerHeight)
A_c = 2 #Square meters, cooling tower cross sectional area (aka Ac)
P_0 = 101325 #Pa, atmospheric pressure (aka pAir)

dx = 0.2 #Distance step

#Boundary Conditions
T_a_0 = 295.15 #Kelvin, atmopheric air temp AKA T_a(X=0) AKA TaIn
W_a_0 = 0.005 #Humidity Ratio atmospheric AKA W_a(X=0) AKA waIn
T_F_n = 312 #Kelvin, hot water temperature AKA T_F(X=n) aka TfIn
v_dot_W_n = 1 #L/min, hot water volumetric flow rate AKA v_dot_W(X=n)

#property lookups at initial conditions
rhoA = CP.PropsSI("D", "T", T_a_0, "P", P_0, "Air") #kg/m^3, Compute air density at atmospheric conditions
rhoF = CP.PropsSI("D", "T", T_a_0, "P", P_0, "Water")#kg/m^3, Compute water density at atmospheric conditions

mfIn = v_dot_W_n*rhoF/1000/60.0 #hot water mass flow rate, kg/s aka m_dot_W(X=n)

m_dot_A_0 = np.mean([6344,8784]) #kg/h, Rules of thumb states air flow rates between 6344 - 8784 kg/h., using the average expected air flow rate AKA ma
m_dot_A_0 = m_dot_A_0/3600 #kg/s AKA ma

velocity = m_dot_A_0/(rhoA*A_c) #m/s
viscosityAir = CP.PropsSI("V", "T", T_a_0, "P", P_0, "Air") #Pa-s, viscosity of air
k = CP.PropsSI("CONDUCTIVITY", "T", T_a_0, "P", P_0, "Air") #W/m/K, conductivity of air
Pr = CP.PropsSI("Prandtl", "T", T_a_0, "P", P_0, "Air") #Prandtl number
Dab = 0.26e-4 #m^2/s Water air diffusivity atmospheric

#Droplet properties
d_p_0 = 0.05 #m, droplet starting diameter.
v_p_0 = (4/3)*np.pi*(d_p_0/2)**3 #m^3, initial droplet volume
A_p_0 = 4*np.pi*(d_p_0/2)**2 #m^2, initial droplet surface area
massParticle = rhoF*v_p_0 #kg, particle mass
Re = velocity*rhoA*d_p_0/viscosityAir #Reynolds number
mu = CP.HAPropsSI('mu','T',T_a_0,'P',P_0,'W',W_a_0) #pa-s, Dynamic Viscosity
mus = mu
Nu = 2.0 + (0.4*Re**0.5+0.06*Re**(2/3))*Pr**0.4*(mu/mus)**(1/4) #nusselt number
Sc = viscosityAir/Dab #schmidt number
Sh = Nu #Sherwood number (Heat/Mass Transfer Analogy)
h = Nu*k/d_p_0 #heat transfer coefficient
hm = Sh * Dab/d_p_0 #mass transfer coefficient

def towerBalance(errors, inputs):
    m_dot_A_0 = inputs[1]
    TaOut = T_a_0 + 3
    TfOut = T_F_n - 3
    waOut = wSat(TaOut)
    waterDH = mfIn*Hf(T_F_n) - mfIn*0.98*Hf(TfOut)
    airDH = m_dot_A_0*Ha(T_a_0, W_a_0) - m_dot_A_0*Ha(TaOut, waOut)
    errors[1] = waterDH + airDH



def towerSimulation(errors, inputs):
    mf = inputs[1:int(length(inputs)/4)*1]
    Tf = inputs[int(length(inputs)/4)*1+1:int(length(inputs)/4)*2]
    Ta = inputs[int(length(inputs)/4)*2+1:int(length(inputs)/4)*3]
    wa = inputs[int(length(inputs)/4)*3+1:int(length(inputs)/4)*4]

    errors = np.reshape(errors, int(length(inputs)/4), 4)
    for x in range(1,(int(length(inputs)/4)-1)):
        errors[x, 1] = mf[x+1] - mf[x] - hm * A_p_0 * (rhoG(Tf[x+1]) - Rh(Ta[x], wa[x])*rhoG(Ta[x])) * dx*Ac * wa[x]/0.622 /v_p_0
        errors[x, 2] = mf[x+1] - mf[x] + ma * (wa[x] - wa[x+1])
        errors[x, 3] = mf[x+1] * Hf(Tf[x+1]) - mf[x]*Hf(Tf[x]) + (mf[x+1] - mf[x])*Hfg(Tf[x]) + h * A_p_0*(Ta[x] - Tf[x+1])  * dx*Ac * wa[x]/0.622 /v_p_0 - Hfg(Tf[x+1])*hm * A_p_0 * (rhoG(Tf[x+1]) - Rh(Ta[x], wa[x])*rhoG(Ta[x])) * dx*Ac * wa[x]/0.622 /v_p_0
        errors[x, 4] = mf[x+1] * Hf(Tf[x+1]) - mf[x]*Hf(Tf[x]) + (mf[x+1] - mf[x])*Hfg(Tf[x]) + ma * (Ha(Ta[x], wa[x]) - Ha(Ta[x+1], wa[x+1]))
    errors[end, 1] = mf[end] - mfIn
    errors[end, 2] = Tf[end] - T_F_n
    errors[end, 3] = Ta[1] - T_a_0
    errors[end, 4] = wa[1] - W_a_0

    return np.reshape(errors, 1, int(length(inputs))) #need to be reshaped as errors must be returned as a long array

mf = np.linspace(0.0131979, mfIn, int(height/dx))
Tf = np.linspace(T_F_n-3, T_F_n, int(height/dx))
Ta = np.linspace(T_a_0, T_a_0+3, int(height/dx))
wa = np.linspace(W_a_0, 0.00741533, int(height/dx))

inputs0 = np.reshape(np.hstack((mf, Tf, Ta, wa)), int(height/dx)*4)

result = root(towerSimulation, inputs0)

#    m_dot = inputs[0:int(x_n/dx)-1]
#    T = inputs[int(x_n/dx):]

#    errors = np.ones((4, int(x_n/dx)))
#    for distance in range(0, len(m_dot)-1):
#        reynolds = velocity * dx/viscosityAir
#        Sc = viscosityAir/(rhoA*Dab)
#        Sh = 0.43*reynolds**0.58*Sc**0.4
#        hm = Dab*Sh/l
#        rhoFs = CP.PropsSI("D", "T", T[int(distance)], "Q", 1, "Water")

        #State Calculations
    #    errors[0][int(distance)] = (m_dot[int(distance)+1] - m_dot[int(distance)])/dx - hm*A_p_0*(rhoFs - rhoFinf)*A_c/0.622/velocity
    #    errors[1][int(distance)] = (m_dot[int(distance)+1] - m_dot[int(distance)])/dx - m_a*W_a_0 #unsure if this is implemented correctly
    #    errors[2][int(distance)] = ((m_dot[int(distance)+1] - m_dot[int(distance)])/dx)*c_p*(T_a[int(distance)+1] - T_a[int(distance)]) - hm*A_p_0*(T_a_0-T_f[int(distance)])*A_c/0.622/v_p_0
    #    errors[3][int(distance)] = m_a*c_p*(T_a[int(distance)+1] - T_a[int(distance)])-((m_dot[int(distance)+1] - m_dot[int(distance)])/dx)*c_p*(T_a[int(distance)+1]

        #Boundary Conditions
    #    errors[0][-1] = m_dot[-1] - m_dot_W_n
    #    errors[1][-1] = m_a[-1] - (v_total - (W_a_0*(v_total))
    #    errors[2][-1] = T_f[-1] - t_outlet_water
    #    errors[3][-1] = T_a[-1] - T_a_0
