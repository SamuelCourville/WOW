from eq36 import *

class aqStep:
    def callAqRockEquil(rockDict, waterDict, press, temp, WR):
        #calcWR() 
        heat=0
        if 'H2O2' in waterDict: 
            return runEq36(temp, press, waterDict, rockDict,WR)
        return rockDict,waterDict,heat

    def callAqEquil(waterDict,pH2,pH,press,temp,name,ind):
        heat=0
        if checkEQfileExist(name,ind):
            pH, pH2, precip, aqSpec, aqGas, heat = runEq6(pH2,pH,temp,press,waterDict,name,str(ind))
        else:
            innitEQfiles(name,str(ind)) 
            pH, pH2, precip, aqSpec, aqGas, heat = runEq36(pH2,pH,temp,press,waterDict,name,str(ind))
        return pH, pH2, precip, aqSpec, aqGas, heat

def innitEQfiles(modelName,r):
    copyMasters(modelName,r)

def runEq6(pH2,pH,Temperature,Pressure,reacs,name,ind):
    ni=name+str(ind)
    mH2  = 10**-pH2 
    mCl  = 0.01914 
    ST=Temperature
    r={'Q':0,'S':0,'C':0,'K':0,'N':0,'L':0,'B':0,'F':0,'P':0,'Mn':0}
    w=convert_WttoMol(reacs)
    updateReacts(reacs,name+str(ind))
    #updateRock(r['Q'],r['S'],r['C'],r['K'],r['N'],r['L'],r['B'],r['F'],r['P'],r['Mn'],w['CO2'],0,w['NH3'],0,0,w['CH4'])
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

def runEq36(pH2,pH,Temperature,Pressure,water,name,ind):
    ni=name+str(ind)
    mH2  = 0.03 #10**-pH2 
    mCl  = 0.01914  
    ST=Temperature
    
    w=convert_WttoMol(water)

    rock={'B':0,'Mn':0,'F':0.04,'Q':0.046,'P':0.048,'N':0.0005,'K':0.0001,'L':0.001,'C':0.003,'S':0.02}

    #rock={'Q':0,'S':0,'C':0,'K':0,'N':0,'L':0,'B':0,'F':0,'P':0,'Mn':0}
    #for i in rockFrac:
    #    rock[i]=convert_VoltoMol(i,rockFrac[i])
    
    updateRock(rock['Q'],rock['S'],rock['C'],rock['K'],rock['N'],rock['L'],rock['B'],rock['F'],rock['P'],rock['Mn'],w['CO2'],0,w['NH3'],0,0,w['CH4'],name+str(ind))
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
    pH =extractpH2(Fugs)
    heat=1

    #print(endMin)
    #print(endAq)
    #print(endFugs)
    
    for i in endMin:
        endMin[i]=(10**endMin[i])/55.55
 
    #return aqout, oceanAQs, solidIce, Xi, Aqs, minerals, Fugs
    return pH, pH2, endMin, endAq, endFugs, heat

def convert_WttoMol(w):
    newW={}
    H2O_M=w['H2O']
    sumM=0
    for i in w:
        sumM+=w[i]
    for i in w:
        if not i=="H2O":
            newW[i]=w[i]/sumM*55.55 # FIX!!!!!
    return newW      # yeah, this is all wrong. Check molarity vs molality 
