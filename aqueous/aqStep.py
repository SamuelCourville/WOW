from eq36 import *
import numpy as np

class aqStep:
    def callAqRockEquil(rockDict, waterDict, press, temp, WR, name, ind, M_Cl):
        #calcWR() 
        heat=0
        if waterDict["H"]>0 and waterDict["O"]>0:
            innitEQfiles(name, str(ind))
            pH, pH2, precip, aqSpec, aqAct, aqGas, heat = runEq36_rock(temp, press, waterDict, rockDict, WR, name, str(ind), M_Cl)
            return aqSpec, aqAct, pH, pH2, aqGas, precip
        else:
            return {}, {}, 0, 0, {}, {}

    def callAqEquil(waterDict,press,temp,name,ind, M_Cl):
        heat=0
        #ElstoAdd={}
        #elementList=["H","C","S","N","Mg","Si","Fe","Ca","Na","Al","K","O"]
        #for i in elementList:
        #    ElstoAdd[i]=0
        #    if i in waterDict:
        #        ElstoAdd[i]=waterDict[i]
        #        waterDict[i]=0
        if temp<263:#checkEQfileExist(name,ind): # use FREZchem here for brine
            aqSpec={}
            aqGas={}
            aqAct={}
            pH=np.nan
            pH2=np.nan
            precip={}
            #pH, pH2, precip, aqSpec, aqGas, heat = runEq6(temp,press,waterDict,name,str(ind))
        else:
            innitEQfiles(name,str(ind)) 
            pH, pH2, precip, aqSpec, aqAct, aqGas, heat = runEq36(temp,press,waterDict,name,str(ind),M_Cl)
        return pH, pH2, precip, aqSpec, aqAct, aqGas, heat

    def inv_refactor_water(dict):
        mH2O=dict["H2O"]
        if "H" in dict:
            dict["H"] += 1 / 9 * mH2O
        else:
            dict["H"] = 1 / 9 * mH2O
        if "O" in dict:
            dict["O"] += 8 / 9 * mH2O
        else:
            dict["O"] = 8 / 9 * mH2O
        del dict['H2O']
        return dict

    def refactor_water(AqComp):

        sumM=0
        for i in AqComp:
            sumM+=AqComp[i]

        AqCompNew={}

        if "H" not in AqComp:
            AqCompNew["H"]=0
        if "O" not in AqComp:
            AqCompNew["O"]=0
        for i in AqComp:
            AqCompNew[i]=AqComp[i]
        mFH=2/18
        mFO=16/18
        if AqCompNew["H"]>mFH*(AqCompNew["H"]+AqCompNew["O"]):
            mH = mFH*(AqCompNew["H"]+AqCompNew["O"])
            AqCompNew["H"] = AqCompNew["H"]-mH
            mO = AqCompNew["O"]
            AqCompNew["O"]=0
            if "H2O" in AqCompNew:
                AqCompNew["H2O"]=AqCompNew["H2O"]+mO+mH
            else:
                AqCompNew["H2O"] = mO + mH
        else:
            mO = mFO*(AqCompNew["H"]+AqCompNew["O"])
            AqCompNew["O"] = AqCompNew["O"]-mO
            mH = AqCompNew["H"]
            AqCompNew["H"]=0
            if "H2O" in AqCompNew:
                AqCompNew["H2O"]=AqCompNew["H2O"]+mO+mH
            else:
                AqCompNew["H2O"] = mO + mH

        sumM2 = 0
        for i in AqCompNew:
            sumM2 += AqCompNew[i]

        if abs(sumM2-sumM)>100000:
            print("failed to keep mass constant")
            print(sumM)
            print(sumM2)
            print(AqCompNew)
            print(AqComp)
            print(error)

        return AqCompNew

    def simple_freeze(AqComp, IceCompNew, temp):
        heat=0
        frozen=False
        IceMeltHeat = 334 * 1000  # J/kg
        if temp<calcFreezePoint(AqComp):
            frozen=True
            heat = AqComp["H2O"] * IceMeltHeat
        return heat, frozen, IceCompNew

    # implement Frezchem
    def freeze_complex(cell,AqComp,AqOriginal,IceComp,press,temp,name,ind):
        AqCompNew={}
        frozen=False
        partial=False
        heat = 0
        fp=calcFreezePoint(AqComp)
        newTemp=temp
        if temp<fp:
            dT=fp-temp
            newTemp=fp
            dE=dT*cell.Cp*cell.Mass
            IceMeltHeat = 334 * 1000  # J/kg

            mH2O_frozen=dE/IceMeltHeat

            if AqComp["H2O"]>mH2O_frozen:
                AqCompNew = AqOriginal
                AqCompNew["H"]=AqCompNew["H"]-1/9*mH2O_frozen
                AqCompNew["O"] = AqCompNew["O"] - 8 / 9 * mH2O_frozen
                if "H" in IceComp:
                    IceComp["H"] += 1/9*mH2O_frozen
                else:
                    IceComp["H"] = 1/9*mH2O_frozen
                if "O" in IceComp:
                    IceComp["O"] += 8/9*mH2O_frozen
                else:
                    IceComp["O"] = 8/9*mH2O_frozen
                partial=True
            else:
                #dT = dE / cell.Cp / cell.Mass
                #newTemp=fp-dT
                iceH=IceComp["H"]
                iceO=IceComp["O"]
                frozen=True
                if "H2O" in IceComp:
                    IceComp["H2O"] += AqComp["H2O"]
                else:
                    IceComp["H2O"] = AqComp["H2O"]
                testAdd=0
                for i in AqComp:
                    if not i=="H2O":
                        IceComp["H2O"]+=AqComp[i]
                        testAdd+=AqComp[i]
                        AqCompNew[i]=AqComp[i]

                heat = AqComp["H2O"] * IceMeltHeat
                if "H2O" in IceComp:
                    IceComp = aqStep.inv_refactor_water(IceComp)
                if "H2O" in AqCompNew:
                    AqCompNew = aqStep.inv_refactor_water(AqCompNew)

                AqCompNew["H"]=-(IceComp["H"]-iceH-AqOriginal['H'])
                AqCompNew["O"]=-(IceComp["O"]-iceO-AqOriginal['O'])

        return heat, frozen, partial, IceComp, AqCompNew, newTemp



    def freeze(AqComp,AqOriginal,press,temp,name,ind):
        IceComp={}
        AqCompNew={}
        frozen=False
        heat = 0
        if temp<calcFreezePoint(AqComp):
            frozen=True
            IceComp["H2O"]=AqComp["H2O"]
            testAdd=0
            for i in AqComp:
                if not i=="H2O":
                    IceComp["H2O"]+=AqComp[i]
                    testAdd+=AqComp[i]
                    AqCompNew[i]=AqComp[i]
            #print("H then O")
            #print(1/9*testAdd)
            #print(8/9*testAdd)
            IceMeltHeat = 334 * 1000  # J/kg
            heat = IceComp["H2O"] * IceMeltHeat
            if "H2O" in IceComp:
                IceComp = aqStep.inv_refactor_water(IceComp)
            if "H2O" in AqCompNew:
                AqCompNew = aqStep.inv_refactor_water(AqCompNew)

            AqCompNew["H"]=-(IceComp["H"]-AqOriginal['H'])
            AqCompNew["O"]=-(IceComp["O"]-AqOriginal['O'])

        return heat, frozen, IceComp, AqCompNew


def calcFreezePoint(AqComp):
    mols=convert_WttoMol(AqComp)
    tot=0
    for i in mols:
        if not i=="H2O":
            tot+=mols[i]
    Kf=1.86
    dT=tot*Kf
    fT=273.15-dT
    return fT


def innitEQfiles(modelName,r):
    copyMasters(modelName,r)



# Check weight units of input-output
# Correct input of Cl and pH2
# Water mixing
# redundant water mass? In dictionary and mWater? 

# Need to extract heat if its worth it
def runEq6(Temperature,Pressure,reacs,name,ind):
    ni=name+str(ind)
    ST=Temperature
    w=convert_WttoMol(reacs)
    updateSpecialReactant(w,ni)
    #updateReacts(reacs,name+str(ind)) # DELETE
    executeEQ6(Temperature,Pressure,name,str(ind))
    N,Xi = numberOfSteps(ni)
    minerals = extractProduct(N,ni)
    Aqs,Acts = extractAq(N,ni)
    Fugs = extractFug(N,ni)
    endAq = extractEnd(Aqs)
    endActs = extractEnd(Acts)
    endMin = extractEnd(minerals)
    endFugs=extractEnd(Fugs)
    pH =extractpH(ni)
    pH2=extractpH2(Fugs)
    heat=1
    return pH, pH2, endMin, endAq, endActs, endFugs, heat

def extractpH2(dic):
    if 'H2(aq)' in dic:
    	return -np.log10(dic['H2(aq)'])
    return 0

def runEq36(Temperature,Pressure,water,name,ind,M_Cl):
    ni=name+str(ind)
    mH2  = 0.0000000000 #10**-pH2
    mCl  = M_Cl
    ST=Temperature
    
    w=convert_WttoMol(water)

    updateSpecialReactant(w,ni)
    updateAQ(mH2,mCl,name+str(ind))
    
    runModel(Temperature,Pressure,name,str(ind))
    N,Xi = numberOfSteps(ni)
    minerals = extractProduct(N,ni)
    Aqs,Acts = extractAq(N,ni)
    Fugs = extractFug(N,ni)
    endAq = extractEnd(Aqs)
    endActs = extractEnd(Acts)
    endMin = extractEnd(minerals)
    endFugs=extractEnd(Fugs)
    pH =extractpH(ni)
    pH2 =extractpH2(Fugs)
    heat=1

    for i in endMin:
        endMin[i]=(10**endMin[i])/55.55
 
    return pH, pH2, endMin, endAq, endActs, endFugs, heat


def runEq36_rock(Temperature, Pressure, water, rock, WR, name, ind, M_Cl):
    ni = name + str(ind)
    mH2 = 0.0000000000
    mCl = M_Cl

    refacted_w=aqStep.refactor_water(water)
    Wmass=refacted_w["H2O"]
    Rmass = Wmass / WR

    sum_w=0
    for i in water:
        sum_w+=water[i]

    sum_r=0
    for i in rock:
        sum_r+=rock[i]

    react_wts={}
    for i in refacted_w:
        if not i=="H2O":
            react_wts[i]=refacted_w[i]/sum_w
    if sum_r>0:
        for i in rock:
            if i in react_wts:
                react_wts[i]+=rock[i]/sum_r/WR
            else:
                react_wts[i] = rock[i]/sum_r/WR

    w = convert_norm_wt_to_Mol(react_wts)

    updateSpecialReactant(w, ni)
    updateAQ(mH2, mCl, name + str(ind))
    runModel(Temperature, Pressure, name, str(ind))
    N, Xi = numberOfSteps(ni)
    minerals = extractProduct(N, ni)
    Aqs, Acts = extractAq(N, ni)
    Fugs = extractFug(N, ni)
    endAq = extractEnd(Aqs)
    endActs = extractEnd(Acts)
    endMin = extractEnd(minerals)
    endFugs = extractEnd(Fugs)
    pH = extractpH(ni)
    pH2 = extractpH2(Fugs)
    heat = 1

    for i in endMin:
        endMin[i] = (10 ** endMin[i]) / 55.55

    return pH, pH2, endMin, endAq, endActs, endFugs, heat

def convert_WttoMol(w):
    #Masses={"H":1,"C":12,"S":32,"N":28,"Mg":24.3,"Si":28.1,"Fe":55.8,"Ca":48.1,"Na":23,"Al":27,"K":39.1,"O":16,"H2O":18}
    Masses={"H":1,"C":12.01,"Mg":24.31,"Si":28.09,"Fe":55.85,"Ca":40.08,"Na":22.99,"Al":26.98,"K":39.10,"O":16,"S":32.07,"N":14.01, "H2O":18}    
    newW={}
    sumM=0
    for i in w:
        sumM+=w[i]
    for i in w:
        newW[i]=(w[i]/sumM)/(Masses[i]/1000) # Check!
    return newW      # yeah, this is wrong. Check molarity vs molality

def convert_norm_wt_to_Mol(w):
    #Masses={"H":1,"C":12,"S":32,"N":28,"Mg":24.3,"Si":28.1,"Fe":55.8,"Ca":48.1,"Na":23,"Al":27,"K":39.1,"O":16,"H2O":18}
    Masses={"H":1,"C":12.01,"Mg":24.31,"Si":28.09,"Fe":55.85,"Ca":40.08,"Na":22.99,"Al":26.98,"K":39.10,"O":16,"S":32.07,"N":14.01}
    newW={}
    for i in w:
        if i in Masses:
            newW[i]=w[i]/Masses[i]*55.55 #(w[i]/sumM)/(/1000) # Check!
    return newW      # yeah, this is wrong. Check molarity vs molality
