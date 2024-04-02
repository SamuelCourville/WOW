from eq36 import *
import numpy as np

class aqStep:
    def callAqRockEquil(rockDict, waterDict, press, temp, WR):
        #calcWR() 
        heat=0
        if 'H2O2' in waterDict: 
            return runEq36(temp, press, waterDict, rockDict,WR)
        return rockDict,waterDict,heat

    def callAqEquil(waterDict,press,temp,name,ind):
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
            pH=np.nan
            pH2=np.nan
            precip={}
            #pH, pH2, precip, aqSpec, aqGas, heat = runEq6(temp,press,waterDict,name,str(ind))
        else:
            innitEQfiles(name,str(ind)) 
            pH, pH2, precip, aqSpec, aqGas, heat = runEq36(temp,press,waterDict,name,str(ind))
        return pH, pH2, precip, aqSpec, aqGas, heat

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

        if abs(sumM2-sumM)>10000:
            print("failed to keep mass constant")
            print(error)

        return AqCompNew

    # implement Frezchem
    def freeze(AqComp,press,temp,name,ind):
        IceComp={}
        AqCompNew={}
        frozen=False
        heat = 0
        if temp<calcFreezePoint(AqComp):
            frozen=True
            IceComp["H2O"]=AqComp["H2O"]
            for i in AqComp:
                if not i=="H2O":
                    IceComp["H2O"]+=AqComp[i]
                    AqCompNew[i]=AqComp[i]
            IceMeltHeat = 334 * 1000  # J/kg
            heat = IceComp["H2O"] * IceMeltHeat
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
    Aqs = extractAq(N,ni)
    Fugs = extractFug(N,ni)
    endAq = extractEnd(Aqs)
    endMin = extractEnd(minerals)
    endFugs=extractEnd(Fugs)
    pH =extractpH(ni)
    pH2=extractpH2(Fugs)
    heat=1
    return pH, pH2, endMin, endAq, endFugs, heat

def extractpH2(dic):
    if 'H2(aq)' in dic:
    	return -np.log10(dic['H2(aq)'])
    return 0

def runEq36(Temperature,Pressure,water,name,ind):
    ni=name+str(ind)
    mH2  = 0.0000000000 #10**-pH2
    mCl  = 0.0001914
    ST=Temperature
    
    w=convert_WttoMol(water)

    updateSpecialReactant(w,ni)
    updateAQ(mH2,mCl,name+str(ind))
    
    runModel(Temperature,Pressure,name,str(ind))
    N,Xi = numberOfSteps(ni)
    minerals = extractProduct(N,ni)
    Aqs = extractAq(N,ni)
    Fugs = extractFug(N,ni)
    endAq = extractEnd(Aqs)
    endMin = extractEnd(minerals)
    endFugs=extractEnd(Fugs)
    pH =extractpH(ni)
    pH2 =extractpH2(Fugs)
    heat=1

    for i in endMin:
        endMin[i]=(10**endMin[i])/55.55
 
    return pH, pH2, endMin, endAq, endFugs, heat

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
