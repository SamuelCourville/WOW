import csv


class LookupProps:
    def calcThermalCond(P,T,RockPhases,RockPhaseDat,IceComp,AqComp,M,rockComp):
        valsR_pre=[]
        wtsR_pre=[]
        if len(RockPhases) == 0:
            for h in rockComp:
                key = LookupProps.crossRef('and_rs')
                wtsR_pre.append(rockComp[h]/M)
                C = LookupProps.getTcondCoeffs(key)
                valsR_pre.append(LookupProps.TCondFunc(C,P,T))
        valsR=[]
        wtsR=[]
        for i,j in enumerate(RockPhases):
            if j not in ['Bulk_rs']:
                key = LookupProps.crossRef(j)
                C = LookupProps.getTcondCoeffs(key)
                wtsR.append(RockPhaseDat[i]['wt%']/100)
                valsR.append(LookupProps.TCondFunc(C,P,T))
        valsI=[]
        wtsI=[]
        for i,k in enumerate(IceComp):
            wtsI.append(IceComp[k]/M)
            C = LookupProps.getTcondCoeffs("H2O") #Assumes everything is H2O! k)
            valsI.append(LookupProps.TCondFunc(C,P,T))
        valsA=[]
        wtsA=[]
        for i,l in enumerate(AqComp):
            wtsA.append(AqComp[l]/M)
            C = LookupProps.getTcondCoeffs("Liquid water") #WRONG! Assumes everything is water.
            valsA.append(LookupProps.TCondFunc(C,P,T))
        vals=valsR+valsI+valsA+valsR_pre
        wts=wtsR+wtsI+wtsA+wtsR_pre
        #print('k:')
        #print(vals)
        #print(wts) 
        tcond = LookupProps.avgProps(vals,wts)
        return tcond


    def calcHeatCap(P,T,RockPhases,RockPhaseDat,IceComp,AqComp,M,rockComp):
        valsR_pre=[]
        wtsR_pre=[]
        if len(RockPhases) == 0:
            for h in rockComp:
                key = LookupProps.crossRef('and_rs')
                wtsR_pre.append(rockComp[h]/M)
                C = LookupProps.getCpCoeffs(key)
                valsR_pre.append(LookupProps.HeatCapFunc(C,P,T))
        valsR=[]
        wtsR=[]
        for i,j in enumerate(RockPhases):
            if j in ['Bulk_rs']:
                #print(j)
                #key = LookupProps.crossRef(j)
                wtsR.append(RockPhaseDat[i]['wt%']/100)
                valsR.append(RockPhaseDat[i]['Heat Capacity (J/K/kg)'])
        valsI=[]
        wtsI=[]
        for i,k in enumerate(IceComp):
            wtsI.append(IceComp[k]/M)
            C = LookupProps.getCpCoeffs("H2O") #Assumes water k)
            valsI.append(LookupProps.HeatCapFunc(C,P,T))
        valsA=[]
        wtsA=[]
        for i,l in enumerate(AqComp):
            wtsA.append(AqComp[l]/M)
            C = LookupProps.getCpCoeffs("Liquid water") # Wrong! assumes everything is water
            valsA.append(LookupProps.HeatCapFunc(C,P,T))
        vals=valsR+valsI+valsA+valsR_pre
        wts=wtsR+wtsI+wtsA+wtsR_pre
        cp = LookupProps.avgProps(vals,wts)
        #print('Cp:')
        #print(vals)
        #print(wts) 
        #print(cp) 
        return cp


    def calcDens(P,T,RockPhases,RockPhaseDat,IceComp,AqComp,M,rockComp):
        valsR_pre=[]
        wtsR_pre=[]
        if len(RockPhases) == 0:
            for h in rockComp:
                key = LookupProps.crossRef('and_rs')
                wtsR_pre.append(rockComp[h]/M)
                C = LookupProps.getRhoCoeffs(key)
                valsR_pre.append(LookupProps.rhoFunc(C,P,T))
        valsR=[]
        wtsR=[]
        for i,j in enumerate(RockPhases):
            if j not in ['Bulk_rs']:
                key = LookupProps.crossRef(j)
                wtsR.append(RockPhaseDat[i]['wt%']/100)
                valsR.append(RockPhaseDat[i]['Density(kg/m3)'])
        valsI=[]
        wtsI=[]
        for i,k in enumerate(IceComp):
            wtsI.append(IceComp[k]/M)
            C = LookupProps.getRhoCoeffs("H2O") # assumes water
            valsI.append(LookupProps.rhoFunc(C,P,T))
        valsA=[]
        wtsA=[]
        for i,l in enumerate(AqComp):
            wtsA.append(AqComp[l]/M)
            C = LookupProps.getRhoCoeffs("Liquid water") # Wrong! assumes everything is water
            valsA.append(LookupProps.rhoFunc(C,P,T))
        vals=valsR+valsI+valsA+valsR_pre
        wts=wtsR+wtsI+wtsA+wtsR_pre
        rho = LookupProps.avgProps(vals,wts)
        #print('rho:')
        #print(vals)
        #print(wts)
        #print(rho) 
        return rho





    def HeatCapFunc(C, P, T):
        val = C[0]+C[1]*T+C[2]*T**2+C[3]*T**3+C[4]*T**-2
        return val

    def TCondFunc(C, P, T):
        val = C[0]+C[1]/T+C[2]*T
        return val

    def rhoFunc(C, P, T):
        val = C[0]
        return val



    
    def avgProps(vals,wts):
        avgVal=0
        for i in range(0,len(vals)):
            avgVal=avgVal+wts[i]*vals[i]
        return avgVal

    def crossRef(spec):
        spec = spec[:-3] # removes '_rs' from perplex name 
        cross_dict=LookupProps.read_2col_csv_to_dict("/Users/samuelcourville/Documents/JPL/combinedModel/tables/perplex_lookup.csv")
        if spec in cross_dict:
            key = cross_dict[spec]
        else:
            print("The mineral "+spec+" is not in my tables. Using olivine instead.")
            key = "Sillicates"
        return key

    def getTcondCoeffs(key):
        coDict=LookupProps.read_Ncol_csv_to_dict("/Users/samuelcourville/Documents/JPL/combinedModel/tables/k.csv")
        return coDict[key]

    def getCpCoeffs(key):
        coDict=LookupProps.read_Ncol_csv_to_dict("/Users/samuelcourville/Documents/JPL/combinedModel/tables/Cp.csv")
        return coDict[key]

    def getRhoCoeffs(key):
        coDict=LookupProps.read_Ncol_csv_to_dict("/Users/samuelcourville/Documents/JPL/combinedModel/tables/rho.csv")
        return coDict[key]

    def read_2col_csv_to_dict(file_path):
        result_dict = {}
        with open(file_path, 'r') as file:
            csv_reader = csv.reader(file)
            for row in csv_reader:
                if len(row) == 2:
                    key, value = row
                    result_dict[key] = value
        return result_dict

    def read_Ncol_csv_to_dict(file_path):
        result_dict = {}
        first_row = True  # Flag to track the first row
        with open(file_path, 'r') as file:
            csv_reader = csv.reader(file)
            for row in csv_reader:
                if first_row:
                    first_row = False  # Skip the first row
                    continue
                if row:  # Ensure the row is not empty
                    key = row[0]
                    values = [float(i) for i in row[1:]]
                    result_dict[key] = values    
        return result_dict


