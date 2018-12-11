#Cooling Tower Model
#Shreyas Sudhakar
#2018

using ExcelReaders
using DifferentialEquations
using NLsolve
using Plots
using PyCall
using CSV
using Statistics
@pyimport CoolProp.CoolProp as CP

energyUse = []
function getNetPresentValue()
    global energyuse
    annualTowerResults = runAnnualSimulation(7) # runAnnualSimulation(365)
    annualDemand = runAnnualDemand(7)
    energyuse = calculateEnergyUse(annualTowerResults)
    annualCosts = calculateAnnualCosts(energyuse)
    capitalCosts = calculateCapitalCosts(mfIn, hxD, hxL)
    npv = calculateNetPresentValue(annualCosts, capitalCosts)
    return npv
end

function waterPumpDeltaP(mdotw, pD, pR, pL)
    velocity = mdotw/(rhoF*pi*pD^2/4)
    viscosityF = 0.001 #N*s/m^2, water viscosity
    Re = velocity*pD/viscosityF

    if Re < 0.0001
        return 0
    else

    A=(2.457*log10(1/((7/Re)^0.9+0.27*pR/pD)))^16
    B=(37530/Re)^16
    f=8*((8/Re)^12+1/(A+B)^(3/2))^(1/12)
    f*pL*velocity^2/(pD*2*gc)
    end
end

function coolingTowerPerformance(makeupwater, hxInT, hxOutT, airflowrate)
    #(* Fan and Pump (makeup, main) Power (kWh) *)
    days = 1
    result = [
        fanDeltaP*airflowrate *days*24/etaF/1000,
        makeupwater /rhoF*9.81*waterPumpDeltaP(makeupwater,pD,pR,pL) * days*24/etaP/1000,
        mfIn /rhoF*9.81*waterPumpDeltaP(mfIn,pD,pR,pL)*days*24/etaP/1000,
    # (* Heat Exchanger Load (kWh) *)
        Ch*(hxOutT-hxInT)*days*24 /1000]
    return result
end

function calculateEnergyUse(dailyEnergyUse)
    totalEnergyUse = []
    (makeup, hxColdInletT, hxColdOutletT, airFlow) = dailyEnergyUse
    for day in 1:length(makeup)
        push!(totalEnergyUse, coolingTowerPerformance(
            makeup[day],
            hxColdInletT[day],
            hxColdOutletT[day],
            airFlow[day]
        ))
    end
    totalEnergyUse
end

#Economic Parameters
electricityCost = 0.054 #$ per kWh, electricity cost
heatLoadValue = 0.025 #$ per kWh, heat load value
watercost = 0.23 #$ per m^3, cost of makeup water

discount = 0.1 #discount rate

function calculateAnnualCosts(dailyEnergyUse)
    fanCosts = 0
    makeuppumpCosts = 0
    waterpumpCosts = 0
    heatLoadCosts = 0
    for day in 1:length(dailyEnergyUse)
        (fanPower, makeuppumpPower, waterpumpPower, heatLoad) = dailyEnergyUse[day]

        fanCosts += -fanPower*electricityCost
        makeuppumpCosts += -makeuppumpPower*electricityCost
        waterpumpCosts += -waterpumpPower*electricityCost
        heatLoadCosts += heatLoad*heatLoadValue
    end
    return [fanCosts, makeuppumpCosts, waterpumpCosts, heatLoadCosts]
end

function calculateNetPresentValue(annualCosts, capitalCost)
    #(* Economic Parameters *)
    projectLife=25
    inflation=0.04 #inflation rate
    fuelInflation=0.013 #Fuel inflation rate

    (fanCosts, makeuppumpCosts, waterpumpCosts, heatLoadCosts) = annualCosts
    npv = 0
    for year in 1:projectLife
        npv += fanCosts/(1+inflation)^year +
            (makeuppumpCosts + waterpumpCosts)/(1+inflation)^year +
            heatLoadCosts/(1+fuelInflation)^year
    end
    npv -= sum(capitalCost)
    return npv
end

capexTower(mfIn=1) = 120000*(mfIn/570)^0.65
capexHX(diameter,length) =500*(3.1415*diameter*length)^0.65
function calculateCapitalCosts(mfIn, hxD, hxL)
    towerCost = capexTower(mfIn*60)
    hxCost = 500*(3.1416*hxL*hxD)^0.65

    return [towerCost, hxCost]
end

function runAnnualSimulation(days=7)
    global dailyData
    xls = openxl("seattledailydata.xlsx") #open excel steam table file
    dailyData = readxlsheet(xls, "seattledailydata") #read in properties as an array

    day = 0
    makeup = zeros(days)
    airFlow = zeros(days)
    hxColdInletT = zeros(days)
    hxColdOutletT = zeros(days)
    for date in 1:(Int(length(dailyData)/4-1))
        #global day, makeup, hxColdInletT, hxColdOutletT
        day = day + 1

        rh = (dailyData[day+1, 2])/100 #relative humidity, percet
        TaIn = dailyData[day+1, 3] + 273.15 #Kelvin, atmopheric air temp AKA T_a(X=0) AKA TaIn
        waIn = omega(TaIn, rh)
        #println(day, "; ", TaIn, "; ", rh)

        mf = LinRange(mfIn*0.999, mfIn, Int(towerHeight/dx))
        Tf = LinRange(TfIn, TfIn, Int(towerHeight/dx))
        Ta = LinRange(TaIn, TaIn, Int(towerHeight/dx))
        wa = LinRange(waIn, waIn*1.001, Int(towerHeight/dx))
        inputs0 = reshape(hcat(mf, Tf, Ta, wa), 1, Int(towerHeight/dx)*4)

        try
            result = nlsolve((errors, inputs) -> towerSimulation(errors, inputs, TaIn, waIn), inputs0,  method=:trust_region, iterations=1000, xtol=0.01)
            solution = result.zero

            mf = solution[1:Int(length(solution)/4)*1]
            Tf = solution[Int(length(solution)/4)*1+1:Int(length(solution)/4)*2]
            Ta = solution[Int(length(solution)/4)*2+1:Int(length(solution)/4)*3]
            wa = solution[Int(length(solution)/4)*3+1:Int(length(solution)/4)*4]

            #println("Day: ", day, " TfOut: ", Tf[1], " TaOut: ", Ta[end])

            airFlow[day] = ma/rhoA + ma*wa[end]/rhoF
            makeup[day]=(mf[end] - mf[1])
            hxColdInletT[day] = (makeup[end]*TfmakeupIn + mf[1]*Tf[1])/mf[end]

            hxColdOutletT[day] =  (Cc/(hxU*hxA)*hxColdInletT[day]+hxHotT)/(1.0+Cc/(hxU*hxA))


        catch errorMessage
            println("Day: ", day, " did not converge. ", errorMessage)
        end

    end
    return [makeup, hxColdInletT, hxColdOutletT, airFlow]
end

function runAnnualDemand(days=7)
    day = 0
    demand = zeros(days)
    for date in 1:(Int(length(dailyData)/4-1))
        #global day, makeup, hxColdInletT, hxColdOutletT
        day = day + 1
        TaIn = dailyData[day+1, 3] + 273.15 #Kelvin, atmopheric air temp AKA T_a(X=0) AKA TaIn
        GHI = dailyData[day+1, 4] #W/m2, Global Horizontal Irradiance
        demand[day] = demand(TaIn,GHI)
    end
    return demand
end

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

function omega(T, rh)
    CP.HAPropsSI("W", "T", T, "P", 101325, "R", rh) #compute humidity ratio for any relative humidity
end

#properties at initial conditions
rhoF = 998.9 #kg/m^3, Compute water density at atmospheric conditions, using constant to save computation
rhoA = 1.103 #kg/m^3, air density at atmospheric conditions, using constant to save computation
P_0 = 101325 #Pa, atmospheric pressure (aka pAir)


#Tower Parameters
towerHeight = 1 #Meters, cooling tower Height
A_c = 2 #Square meters, cooling tower cross sectional area (aka Ac)
dx = 0.1 #Distance step
ma = 150000 #kg/h, THIS IS 18.0 IN PROFS CODE
ma = ma/3600 #kg/s AKA ma
TaIn = 295.15 #Kelvin, atmopheric air temp AKA T_a(X=0) #THIS IS A CONSTANT IN PROFS CODE, fixing as constant to get code to run
waIn = 0.005 #Humidity Ratio atmospheric AKA W_a(X=0) AKA waIn #THIS IS A CONSTANT IN PROFS CODE, fixing as constant to get code to run
v_dot_W_n = 1 #L/min*m^2, hot water volumetric flow rate AKA v_dot_W(X=n) #THIS IS 0.8 IN PROFS CODE
mfIn = v_dot_W_n*rhoF/1000/60.0 #hot water mass flow rate, kg/s*m^2 aka m_dot_W(X=n)
TfIn = 312 #Kelvin, hot water temperature AKA T_F(X=n) aka TfIn

#Heat Exchanger
hxD = 0.025 #m, diameter of heat exchanger tubes
hxL = 20 #meters, length of water pipe in tower
hxA = 2*(pi*hxD/2)*hxL #Surface area of heat exchanger tubes
hxU = 1200 #W/m^2/K, heat transfer coefficient of hex

hxHotT = 323.5 #K, heat exchanger shell side inlet water temperature
vfshellin = 0.75 #L/min, heat exchanger shell side inlet water flow rate
hxHotFlow = vfshellin* rhoF/1000/60 #kg/s, heat exchanger hot water mass flow rate

TfmakeupIn = 288.15 #Kelvin, makeup water temerature

Cpw = CP.PropsSI("C", "T", hxHotT, "P", P_0, "Water")#J/kg/K, Compute CP at atmospheric conditions #THIS IS A CONSTANT IN PROFS CODE
Ch = hxHotFlow*Cpw
Cc = mfIn*Cpw

# Pipe Parameters
gc = 1 #kg * m/N/s^2, acceleration due to gravity
pL = 200 #m, length of pipe
pR = 0.0001 #pipe roughness
pD = 0.1 #m, pipe diameter
etaP = 0.72 #Efficiency of water pump

#Fan
fanDeltaP = 500 #Pa, pressure differential fan needs to overcome
etaF = 0.7 #Efficiency of fan

velocity = ma/(rhoA*A_c) #m/s
viscosityAir = CP.PropsSI("V", "T", TaIn, "P", P_0, "Air") #Pa-s, viscosity of air
k = 0.0279 #W/m/K, conductivity of air
Pr = 0.707 #Prandtl number
Dab = 0.26e-4 #m^2/s Water air diffusivity atmospheric

#Droplet properties
dropletD = 0.05 #m, droplet starting diameter.
volumeParticle = (4/3)*pi*(dropletD/2)^3 #m^3, initial droplet volume
areaParticle = 4*pi*(dropletD/2)^2 #m^2, initial droplet surface area
massParticle = rhoF*volumeParticle #kg, particle mass
Re = velocity*rhoA*dropletD/viscosityAir #Reynolds number
mu = CP.HAPropsSI("mu","T",TaIn,"P",P_0,"W",waIn) #pa-s, Dynamic Viscosity
mus = mu
Nu = 2.0 + (0.4*Re^0.5+0.06*Re^(2/3))*Pr^0.4*(mu/mus)^(1/4) #nusselt number
Sc = viscosityAir/Dab #schmidt number
Sh = Nu #Sherwood number (Heat/Mass Transfer Analogy)
h = Nu*k/dropletD #heat transfer coefficient
hm = Sh * Dab/dropletD #mass transfer coefficient

#Energy Balance
# function towerBalance(errors, inputs)
#     ma = inputs[1]
#     TaOut = TaIn + 3
#     TfOut = TfIn - 3
#     waOut = wSat(TaOut)
#     waterDH = mfIn*Hf(TfIn) - mfIn*0.98*Hf(TfOut)
#     airDH = ma*Ha(TaIn, waIn) - ma*Ha(TaOut, waOut)
#     errors[1] = waterDH + airDH
# end

#result = nlsolve(towerBalance, [ma])
#println(result)

#Building Dimensions
lBuilding = 150 #Feet, length of building
wBuilding = 275 #Feet, width of building
hBuilding = 10 #Feet, height of building
aRoof = lBuilding*wBuilding #ft^2, area of roof (assuming flat roof)
aWalls = (2*(lBuilding+wBuilding))*hBuilding #ft^2, area of walls
vBuilding = aRoof*hBuilding #ft^3, volume of building
vTimber = 0.02 #Volume fraction of building walls made of timber
tWall = 1/2 #inches, thickness of timber walls
absorptivity = 0.35 #Absorptivity of timber
hsummer = 0.09375 #W/m*K, Overall heat transfer coefficient in summer
hwinter = 0.0894 #W/m*K, Overall heat transfer coefficient in winter
vGlass = 0.98 #Volume fraction of building walls made of 2 pane glass
airGap = 3/16 #Inch, air space between panes of glass
transmitGlass = 0.86 #Solar transmission through glass

#Design Parameters
nAirChanges = 12 #Target Number of air changes in commercial building
TInside = 295 #Air conditioner set temperature inside the building
kTimber = 0.15 #W/mK, thermakl conductivty of pine wood, see https://www.engineeringtoolbox.com/thermal-conductivity-d_429.html
kFoamConcrete = 0.2 #W/mK, thermal conductivity of Foam (aka Lightweight) concrete, see https://www.engineeringtoolbox.com/thermal-conductivity-d_429.html
kAir = 0.0262 #W/mK, thermal conductivity of air, see https://www.engineeringtoolbox.com/thermal-conductivity-d_429.html
eConcrete = 0.85 #emmissivity coefficient of concrete, see https://www.engineeringtoolbox.com/emissivity-coefficients-d_447.html

Cpa = CP.PropsSI("C", "T", TInside, "P", P_0, "Air")#J/kg/K, Compute CP at air temperature inside building

function season(T)
    if (T > 280) #Setting as arbitrary temperature threshhold for winter weather
        h = hsummer
    end
    h = hwinter
end

function Qcond(T)
    Q_cond_roof = (kFoamConcrete*aRoof*(TInside - T))/tWall
    Q_cond_timber = (kTimber*(aWalls*vTimber)*(TInside - T))/tWall
    Q_cond_windows = (kAir*(aWalls*vGlass)*(TInside - T))/airGap
    Qcond = Q_cond_roof + Q_cond_timber + Q_cond_windows
end

function Qconv(T)
    h = season(T)
    Qconv = h*(aWalls + aRoof)*(TInside - T)
end

function Qrad(T,GHI)
    Qrad = eConcrete*aRoof*GHI*(TInside^4 - T^4) #Assume that entirety of radiative heat tranfer occurs through roof of building
end

function Qevap(T)
    Qevapventilation = rhoA*((nAirChanges*vBuilding)/3600)*Cpa*(TInside - T)
end

function demand(T,GHI)
    Qevap = Qevap(T)
    Qrad = Qrad(T,GHI)
    Qconv = Qconv(T)
    Qcond = Qcond(T)
    demand = Qevap + Qrad + Qconv + Qcond
end



function towerSimulation(errors, inputs, TaIn, waIn)
    mf = inputs[1:Int(length(inputs)/4)*1]
    Tf = inputs[Int(length(inputs)/4)*1+1:Int(length(inputs)/4)*2]
    Ta = inputs[Int(length(inputs)/4)*2+1:Int(length(inputs)/4)*3]
    wa = inputs[Int(length(inputs)/4)*3+1:Int(length(inputs)/4)*4]

    # println(Ta)
    # println(Tf)
    errors = reshape(errors, Int(length(inputs)/4), 4)
    for x in 1:(Int(length(inputs)/4)-1)
        errors[x, 1] = mf[x+1] - mf[x] - hm * areaParticle * (rhoG(Tf[x+1]) - Rh(Ta[x], wa[x])*rhoG(Ta[x])) * dx*A_c * wa[x]/0.622 /volumeParticle
        errors[x, 2] = mf[x+1] - mf[x] + ma * (wa[x] - wa[x+1])
        errors[x, 3] = mf[x+1] * Hf(Tf[x+1]) - mf[x]*Hf(Tf[x]) + h * areaParticle*(Ta[x] - Tf[x+1])  * dx*A_c * wa[x]/0.622 /volumeParticle + Hfg(Tf[x+1])*hm * areaParticle * (rhoG(Tf[x+1]) - Rh(Ta[x], wa[x])*rhoG(Ta[x])) * dx*A_c * wa[x]/0.622 /volumeParticle
        errors[x, 4] = mf[x+1] * Hf(Tf[x+1]) - mf[x]*Hf(Tf[x]) + (mf[x+1] - mf[x])*Hfg(Tf[x]) + ma * (Ha(Ta[x], wa[x]) - Ha(Ta[x+1], wa[x+1]))
    end
    errors[end, 1] = mf[end] - mfIn
    errors[end, 2] = Tf[end] - TfIn
    errors[end, 3] = Ta[1] - TaIn
    errors[end, 4] = wa[1] - waIn
    # display(errors)

    return reshape(errors, 1, Int(length(inputs))) #need to be reshaped as errors must be returned as a long array
end

getNetPresentValue()

energyuseT = hcat(energyuse...)
fanPower = energyuseT[1, :]
makeuppumpPower = energyuseT[2, :]
waterpumpPower = energyuseT[3, :]
heatLoad = energyuseT[4, :]
days = dailyData[2:end, 1]

plot(days, fanPower, xlabel="Date", ylabel = "Power (kW)", title = "Fan Power", label = "Seattle, WA Boeing Field")
plot(days, makeuppumpPower, xlabel="Date", ylabel = "Power (kW)", title = "Makeup Water Pump Power", label = "Seattle, WA Boeing Field")
plot(days, waterpumpPower, xlabel="Date", ylabel = "Power (kW)", title = "Water Pump Power", label = "Seattle, WA Boeing Field")
plot(days, heatLoad, xlabel="Date", ylabel = "Load (kW)", title = "Heat Load", label = "Seattle, WA Boeing Field")
plot(dailyData[2:end, 1], dailyData[2:end, 2], xlabel="Date", ylabel = "Relative Humidity (%)", title = "Average Daily Relative Humidity", label = "Seattle, WA Boeing Field")
plot(dailyData[2:end, 1], dailyData[2:end, 3], xlabel="Date", ylabel = "Temperature (Degrees C)", title = "Average Daily Dry Bulb Temperature", label = "Seattle, WA Boeing Field")
plot(days, demand, xlabel="Date", ylabel = "Power (kW)", title = "Building Power Demand", label = "Seattle, WA Boeing Field")
