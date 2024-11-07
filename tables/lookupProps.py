import csv
import numpy as np

import WOW

##### UPDATE THIS LINE
main_dir=WOW.main_directory

# Other paths
lookup_dir=main_dir+"tables/perplex_lookup.csv"
k_dir=main_dir+"tables/k.csv"
cp_dir=main_dir+"tables/Cp.csv"
rho_dir=main_dir+"tables/rho.csv"


class LookupProps:
    """
        Class to contain functions related to grabbing data from thermal and physical property tables.
    """
    def loadk():
        """
            Load file with thermal conductivity values
        """
        coDict=LookupProps.read_Ncol_csv_to_dict(k_dir)
        return coDict
    def loadrho():
        """
            Load file with density values
        """
        coDict = LookupProps.read_Ncol_csv_to_dict(rho_dir)
        return coDict
    def loadcp():
        """
            load file with heat capacity values
        """
        coDict = LookupProps.read_Ncol_csv_to_dict(cp_dir)
        return coDict
    def loadperp():
        """
            Load file with cross references between Perple_X minerals and materials in other property tables
        """
        cross_dict = LookupProps.read_2col_csv_to_dict(lookup_dir)
        return cross_dict

    def calcThermalCond(ktab,crosstab,P,T,RockPhases,RockPhaseDat,IceComp,AqComp,M,rockComp):
        """
            Calculate thermal conductivity from table coefficients and P,T
        """
        valsR_pre=[]
        wtsR_pre=[]
        if len(RockPhases) == 0:
            for h in rockComp:
                key = LookupProps.crossRef('undiff_rs',crosstab) ### WHAT IS THIS DOING???
                wtsR_pre.append(rockComp[h]/M)
                C = LookupProps.getTcondCoeffs(key,ktab)
                valsR_pre.append(LookupProps.TCondFunc(C,P,T))
        valsR=[]
        wtsR=[]
        for i,j in enumerate(RockPhases):
            if j not in ['Bulk_rs']:
                key = LookupProps.crossRef(j,crosstab)
                C = LookupProps.getTcondCoeffs(key,ktab)
                wtsR.append(RockPhaseDat[i]['wt%']/100)
                valsR.append(LookupProps.TCondFunc(C,P,T))
        valsI=[]
        wtsI=[]
        for i,k in enumerate(IceComp):
            wtsI.append(IceComp[k]/M)
            C = LookupProps.getTcondCoeffs("Ice",ktab) #Assumes everything is H2O! k)
            valsI.append(LookupProps.TCondFunc(C,P,T))
        valsA=[]
        wtsA=[]
        for i,l in enumerate(AqComp):
            wtsA.append(AqComp[l]/M)
            C = LookupProps.getTcondCoeffs("Liquid water",ktab) #WRONG! Assumes everything is water.
            valsA.append(LookupProps.TCondFunc(C,P,T))
        vals=valsR+valsI+valsA+valsR_pre
        wts=wtsR+wtsI+wtsA+wtsR_pre
        tcond = LookupProps.avgProps(vals,wts)
        return tcond


    def calcHeatCap(cptab,crosstab,P,T,RockPhases,RockPhaseDat,IceComp,AqComp,M,rockComp):
        """
            Calculate heat capacity from table coefficients and P,T
        """
        valsR_pre=[]
        wtsR_pre=[]
        if len(RockPhases) == 0:
            for h in rockComp:
                key = LookupProps.crossRef('undiff_rs',crosstab)
                wtsR_pre.append(rockComp[h]/M)
                C = LookupProps.getCpCoeffs(key,cptab)
                valsR_pre.append(LookupProps.HeatCapFunc(C,P,T))
        valsR=[]
        wtsR=[]
        for i,j in enumerate(RockPhases):
            if j in ['Bulk_rs']:
                wtsR.append(RockPhaseDat[i]['wt%']/100)
                valsR.append(RockPhaseDat[i]['Heat Capacity (J/K/kg)'])
        valsI=[]
        wtsI=[]
        for i,k in enumerate(IceComp):
            wtsI.append(IceComp[k]/M)
            C = LookupProps.getCpCoeffs("H2O",cptab) #Assumes water k)
            valsI.append(LookupProps.HeatCapFunc(C,P,T))
        valsA=[]
        wtsA=[]
        for i,l in enumerate(AqComp):
            wtsA.append(AqComp[l]/M)
            C = LookupProps.getCpCoeffs("Liquid water",cptab) # Wrong! assumes everything is water
            valsA.append(LookupProps.HeatCapFunc(C,P,T))
        vals=valsR+valsI+valsA+valsR_pre
        wts=wtsR+wtsI+wtsA+wtsR_pre
        cp = LookupProps.avgProps(vals,wts)
        return cp

    def getRockDens(rhotab, crosstab, P,T,RockPhases,RockPhaseDat,M,rockComp):
        """
            Calculate density from table coefficients and P,T
        """
        valsR_pre=[]
        wtsR_pre=[]
        if len(RockPhases) == 0:
            for h in rockComp:
                key = LookupProps.crossRef('undiff_rs', crosstab)
                wtsR_pre.append(rockComp[h]/M)
                C = LookupProps.getRhoCoeffs(key,rhotab)
                valsR_pre.append(LookupProps.rhoFunc(C,P,T))
        valsR=[]
        wtsR=[]
        for i,j in enumerate(RockPhases):
            if j not in ['Bulk_rs']:
                key = LookupProps.crossRef(j,crosstab)
                wtsR.append(RockPhaseDat[i]['wt%']/100)
                valsR.append(RockPhaseDat[i]['Density(kg/m3)'])
        vals = valsR + valsR_pre
        wts = wtsR + wtsR_pre

        wts=wts/np.sum(wts)

        nE = len(vals)
        vol = 0
        for i in range(nE):
            vol = vol + wts[i] / vals[i]
        if vol==0:
            return 1.0
        rho = 1.0 / vol
        return rho

    def calcDens(rhotab, crosstab, P,T,RockPhases,RockPhaseDat,IceComp,AqComp,M,rockComp):
        """
            Calculate density from table coefficients and P,T
        """
        valsR_pre=[]
        wtsR_pre=[]
        if len(RockPhases) == 0:
            for h in rockComp:
                key = LookupProps.crossRef('undiff_rs', crosstab)
                wtsR_pre.append(rockComp[h]/M)
                C = LookupProps.getRhoCoeffs(key,rhotab)
                valsR_pre.append(LookupProps.rhoFunc(C,P,T))
        valsR=[]
        wtsR=[]
        for i,j in enumerate(RockPhases):
            if j not in ['Bulk_rs']:
                key = LookupProps.crossRef(j,crosstab)
                wtsR.append(RockPhaseDat[i]['wt%']/100)
                valsR.append(RockPhaseDat[i]['Density(kg/m3)'])
        valsI=[]
        wtsI=[]
        for i,k in enumerate(IceComp):
            wtsI.append(IceComp[k]/M)
            C = LookupProps.getRhoCoeffs("H2O",rhotab) # assumes water
            valsI.append(LookupProps.rhoFunc(C,P,T))
        valsA=[]
        wtsA=[]
        for i,l in enumerate(AqComp):
            wtsA.append(AqComp[l]/M)
            C = LookupProps.getRhoCoeffs("Liquid water",rhotab) # Wrong! assumes everything is water
            valsA.append(LookupProps.rhoFunc(C,P,T))
        vals=valsR+valsI+valsA+valsR_pre
        wts=wtsR+wtsI+wtsA+wtsR_pre

        wts=wts/np.sum(wts)

        nE = len(vals)
        vol = 0
        for i in range(nE):
            vol=vol+wts[i]/vals[i]
        rho = 1.0/vol
        #rho = LookupProps.avgProps(vals,wts)
        #print('rho:')
        #print(vals)
        #print(wts)
        #print(rho) 
        return rho





    def HeatCapFunc(C, P, T):
        """
            Function to calculate heat capacity from coefficients
        """
        val = C[0]+C[1]*T+C[2]*T**2+C[3]*T**3+C[4]*T**-2
        return val

    def TCondFunc(C, P, T):
        """
            Calculate thermal conductivity from table coefficients and P,T
        """
        val = C[0]+C[1]/T+C[2]*T
        return val

    def rhoFunc(C, P, T):
        """
            Calculate thermal conductivity from table coefficients and P,T
        """
        val = C[0]
        return val



    
    def avgProps(vals,wts):
        """
            Weighted average of property values
        """
        avgVal=0
        for i in range(0,len(vals)):
            avgVal=avgVal+wts[i]*vals[i]
        return avgVal

    def crossRef(spec,cross_dict):
        """
            Cross reference perple_x minerals with other table entities.
        """
        spec = spec[:-3] # removes '_rs' from perplex name
        if spec in cross_dict:
            key = cross_dict[spec]
        else:
            print("The mineral "+spec+" is not in my tables. Using olivine instead.")
            key = "Sillicates"
        return key

    def getTcondCoeffs(key,coDict):
        """
            Get k coefficients for specific material
        """
        return coDict[key]

    def getCpCoeffs(key,coDict):
        """
            Get Cp coefficients for specific material
        """
        return coDict[key]

    def getRhoCoeffs(key,coDict):
        """
            Get density coefficients for specific material
        """
        return coDict[key]

    def read_2col_csv_to_dict(file_path):
        """
            Read property csv file into dictionary
        """
        result_dict = {}
        with open(file_path, 'r') as file:
            csv_reader = csv.reader(file)
            for row in csv_reader:
                if len(row) == 2:
                    key, value = row
                    result_dict[key] = value
        return result_dict

    def read_Ncol_csv_to_dict(file_path):
        """
            read property csv file into dictionary
        """
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


