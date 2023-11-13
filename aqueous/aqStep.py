from eq36 import *

class aqStep:
    def callAqRockEquil(rockDict, waterDict, press, temp, WR):
        #calcWR() 
        heat=0
        if 'H2O2' in waterDict: 
            return runEq36(temp, press, waterDict, rockDict,WR)
        return rockDict,waterDict,heat


def runEq36(Temperature,Pressure,water,rockFrac, WR):
   
    mH2  = 0.01 
    mCl  = 0.01914 
    ST=Temperature

    rock={'B':0,'Mn':0,'F':0.04,'Q':0.046,'P':0.048,'N':0.0005,'K':0.0001,'L':0.001,'C':0.003,'S':0.02}

    #rock={'Q':0,'S':0,'C':0,'K':0,'N':0,'L':0,'B':0,'F':0,'P':0,'Mn':0}
    #for i in rockFrac:
    #    rock[i]=convert_VoltoMol(i,rockFrac[i])
    
    updateRock(rock['Q'],rock['S'],rock['C'],rock['K'],rock['N'],rock['L'],rock['B'],rock['F'],rock['P'],rock['Mn'],0,0,0,0,0,0)
    updateAQ(mH2,mCl)
    runModel(Temperature,Pressure)
    N,Xi = numberOfSteps()
    minerals = extractProduct(N)
    Aqs = extractAq(N)
    Fugs = extractFug(N)
    endAq = extractEnd(Aqs)
    endMin = extractEnd(minerals)
    endFugs=extractEnd(Fugs)
    pH =extractpH()
    heat=1

    #print(endMin)
    #print(endAq)
    #print(endFugs)
    
    for i in endMin:
        endMin[i]=(10**endMin[i])/55.55
 
    #return aqout, oceanAQs, solidIce, Xi, Aqs, minerals, Fugs
    return endMin, endAq, heat

def convert_VoltoMol(i,rockFrac):
    return 55.55*rockFrac # yeah, this is all wrong. 
