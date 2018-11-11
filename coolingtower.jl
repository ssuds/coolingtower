#Cooling Tower Model
#Shreyas Sudhakar
#2018

using ExcelReaders
using DifferentialEquations
using NLsolve
using Plots
using PyCall
using Statistics
@pyimport CoolProp.CoolProp as CP

xls = openxl("seattledailydata.xlsx") #open excel steam table file
dailyData = readxlsheet(xls, "seattledailydata") #read in properties as an array

plot(dailyData[2:end, 1], dailyData[2:end, 2], xlabel="Date", ylabel = "Relative Humidity (%)", title = "Average Daily Relative Humidity", label = "Seattle, WA Boeing Field")
plot(dailyData[2:end, 1], dailyData[2:end, 3], xlabel="Date", ylabel = "Temperature (Degrees C)", title = "Average Daily Dry Bulb Temperature", label = "Seattle, WA Boeing Field")


function Hf(T)
    if  (T < 273.15 || T > 645 ) #solution moved out of coolprop bounds
        return -1e6
        print("Hf is Broken!!!!")
    end
    CP.PropsSI("H", "T", T, "P", 101325, "Water") #Computing enthalpy of coolant fluid
end

function Hfg(T)
    if  (T < 273.15 || T > 645 ) #solution moved out of coolprop bounds
        return -1e6
        print("Hfg is Broken!!!!")
    end
    CP.PropsSI("H", "T", T, "Q", 1, "Water") - CP.PropsSI("H", "T", T, "Q", 0, "Water") #Computing enthalpy of vaporization
end


function Psat(T)
    if  (T < 273.15) #solution moved out of coolprop bounds
        return -1e6
        print("Psat is Broken!!!!")
    end
    CP.PropsSI("P", "T", T, "Q", 1, "Water") #computing saturation pressure
end

function rhoG(T)
    if  (T < 273.15 || T > 645 ) #solution moved out of coolprop bounds
        return -1e6
        print("rhoG is Broken!!!!")
    end
    CP.PropsSI("D", "T", T, "Q", 1, "Water") #compute density
end

function wSat(T)
    CP.HAPropsSI("W", "T", T, "P", 101325, "R", 1) #compute humidity ratio for saturated air
end

function Rh(T, w)
    if  (T < 273.15 || T > 645 ) #solution moved out of coolprop bounds
        return -1e6
        print("Rh is Broken!!!!")
    end
    CP.HAPropsSI("R", "T", T, "P", 101325, "W", w) #compute relative humidity given temperature and humidity ratio
end

function Ha(T, w)
    if  (T < 273.15 || T > 645 ) #solution moved out of coolprop bounds
        return -1e6
        print("Ha is Broken!!!!")
    end
    CP.HAPropsSI("H", "T", T, "P", 101325, "W", w) #compute enthalpy of air given temperature and humidity ratio
end



#Define known variables
load = 1000 #Watts, cooling load from business building
height = 3 #Meters, cooling tower height (aka towerHeight)
A_c = 20 #Square meters, cooling tower cross sectional area (aka Ac)
P_0 = 101325 #Pa, atmospheric pressure (aka pAir)

dx = 0.2 #Distance step

#Boundary Conditions
TaIn = 295.15 #Kelvin, atmopheric air temp AKA T_a(X=0) AKA TaIn
waIn = 0.005 #Humidity Ratio atmospheric AKA W_a(X=0) AKA waIn
TfIn = 312 #Kelvin, hot water temperature AKA T_F(X=n) aka TfIn
v_dot_W_n = 1 #L/min*m^2, hot water volumetric flow rate AKA v_dot_W(X=n)

#property lookups at initial conditions
rhoA = CP.PropsSI("D", "T", TaIn, "P", P_0, "Air") #kg/m^3, Compute air density at atmospheric conditions
rhoF = CP.PropsSI("D", "T", TaIn, "P", P_0, "Water")#kg/m^3, Compute water density at atmospheric conditions

mfIn = v_dot_W_n*rhoF/1000/60.0 #hot water mass flow rate, kg/s*m^2 aka m_dot_W(X=n)

maIn = 150000 #kg/h,
maIn = maIn/3600 #kg/s AKA ma

velocity = maIn/(rhoA*A_c) #m/s
viscosityAir = CP.PropsSI("V", "T", TaIn, "P", P_0, "Air") #Pa-s, viscosity of air
k = 0.0279 #W/m/K, conductivity of air
Pr = 0.707 #Prandtl number
Dab = 0.26e-4 #m^2/s Water air diffusivity atmospheric

#Droplet properties
d_p_0 = 0.05 #m, droplet starting diameter.
v_p_0 = (4/3)*pi*(d_p_0/2)^3 #m^3, initial droplet volume
A_p_0 = 4*pi*(d_p_0/2)^2 #m^2, initial droplet surface area
massParticle = rhoF*v_p_0 #kg, particle mass
Re = velocity*rhoA*d_p_0/viscosityAir #Reynolds number
mu = CP.HAPropsSI("mu","T",TaIn,"P",P_0,"W",waIn) #pa-s, Dynamic Viscosity
mus = mu
Nu = 2.0 + (0.4*Re^0.5+0.06*Re^(2/3))*Pr^0.4*(mu/mus)^(1/4) #nusselt number
Sc = viscosityAir/Dab #schmidt number
Sh = Nu #Sherwood number (Heat/Mass Transfer Analogy)
h = Nu*k/d_p_0 #heat transfer coefficient
hm = Sh * Dab/d_p_0 #mass transfer coefficient

# function towerBalance(errors, inputs)
#     maIn = inputs[1]
#     TaOut = TaIn + height
#     TfOut = TfIn - height
#     waOut = wSat(TaOut)
#     waterDH = mfIn*Hf(TfIn) - mfIn*0.98*Hf(TfOut)
#     airDH = maIn*Ha(TaIn, waIn) - maIn*Ha(TaOut, waOut)
#     errors[1] = waterDH + airDH
# end


function towerSimulation(errors, inputs)
    mf = inputs[1:Int(length(inputs)/4)*1]
    Tf = inputs[Int(length(inputs)/4)*1+1:Int(length(inputs)/4)*2]
    Ta = inputs[Int(length(inputs)/4)*2+1:Int(length(inputs)/4)*3]
    wa = inputs[Int(length(inputs)/4)*3+1:Int(length(inputs)/4)*4]

    errors = reshape(errors, Int(length(inputs)/4), 4)
    for x in 1:(Int(length(inputs)/4)-1)
        errors[x, 1] = mf[x+1] - mf[x] - hm * A_p_0 * (rhoG(Tf[x+1]) - Rh(Ta[x], wa[x])*rhoG(Ta[x])) * dx*A_c * wa[x]/0.622 /v_p_0
        errors[x, 2] = mf[x+1] - mf[x] + maIn * (wa[x] - wa[x+1])
        errors[x, 3] = mf[x+1] * Hf(Tf[x+1]) - mf[x]*Hf(Tf[x]) + (mf[x+1] - mf[x])*Hfg(Tf[x]) + h * A_p_0*(Ta[x] - Tf[x+1])  * dx*A_c * wa[x]/0.622 /v_p_0 - Hfg(Tf[x+1])*hm * A_p_0 * (rhoG(Tf[x+1]) - Rh(Ta[x], wa[x])*rhoG(Ta[x])) * dx*A_c * wa[x]/0.622 /v_p_0
        errors[x, 4] = mf[x+1] * Hf(Tf[x+1]) - mf[x]*Hf(Tf[x]) + (mf[x+1] - mf[x])*Hfg(Tf[x]) + maIn * (Ha(Ta[x], wa[x]) - Ha(Ta[x+1], wa[x+1]))
    end
    errors[end, 1] = mf[end] - mfIn
    errors[end, 2] = Tf[end] - TfIn
    errors[end, 3] = Ta[1] - TaIn
    errors[end, 4] = wa[1] - waIn

    return reshape(errors, 1, Int(length(inputs))) #need to be reshaped as errors must be returned as a long array
end

#Initial Guesses
mf = LinRange(0.0131979, mfIn, Int(height/dx))
Tf = LinRange(TfIn-3, TfIn, Int(height/dx))
Ta = LinRange(TaIn, TaIn+3, Int(height/dx))
wa = LinRange(waIn, 0.00741533, Int(height/dx))

inputs0 = reshape(hcat(
    mf,
    Tf,
    Ta,
    wa), 1, Int(height/dx)*4)

result = nlsolve(towerSimulation, inputs0, method=:trust_region, iterations=1000, xtol=0.01)

println(result)

solution = result.zero
print(result.zero)


mf = solution[1:Int(length(solution)/4)*1]
Tf = solution[Int(length(solution)/4)*1+1:Int(length(solution)/4)*2]
Ta = solution[Int(length(solution)/4)*2+1:Int(length(solution)/4)*3]
wa = solution[Int(length(solution)/4)*3+1:Int(length(solution)/4)*4]


plot(mf, ylabel = "Fluid mass, kg", xlabel = "Height, m")
plot(Tf, label="Tf", ylabel = "Temperature, K", xlabel = "Height, m")
plot!(Ta, label="Ta")
plot(wa, ylabel = "Humidity Ratio", xlabel = "Height, m")

println("Solutions")
println(mf)
println(Tf)
println(Ta)
println(wa)

println("Heat Capacity")
println(mf[1]*Hf(Tf[1]) - mf[end]*Hf(Tf[end]))


percent = 1+ mf[Int(height/dx)]-mf[1]/mf[Int(height/dx)]
TfOut = Tf[1]
waOut = wa[Int(height/dx)]
TaOut = Ta[Int(height/dx)]
