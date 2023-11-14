import numpy as np

class Decay:
    ### CONSTANTS
    YR = 3.1536E7   #(seconds)
    LN2 = np.log(2)
    SSAGE=4.571*10**9*YR

    # Heat Production W/kg
    H0U238  = 94.65E-6
    H0U235  = 568.7E-6
    H0TH    = 26.38E-6
    H0K     = 29.17E-6
    H0FE    = 6.7E-2
    H0AL    = 3.57E-1
    H0MN    = 5.E-3
    H0CA    = 1.53E-1

    # Half-Lives
    LAU238  = 4.47E9*YR
    LAU235  = 7.0381E8*YR
    LATH    = 14.01E9*YR
    LAK     = 1.277E9*YR
    LAFE    = 1.5E6*YR
    LAAL    = 7.16E5*YR
    LAMN    = 3.74E6*YR
    LACA    = 0.103E5*YR

    # Put constants in dictionary for easy access
    LAs = dict({'U238': LAU238, 
                  'U235':LAU235,
                  'Th232':LATH,
                  'K40':LAK,
                  'Fe60':LAFE,
                  'Al26':LAAL,
                  'Mn53':LAMN,
                  'Ca41':LACA})
    H0s = dict({'U238': H0U238, 
                  'U235':H0U235,
                  'Th232':H0TH,
                  'K40':H0K,
                  'Fe60':H0FE,
                  'Al26':H0AL,
                  'Mn53':H0MN,
                  'Ca41':H0CA})

    def calcOriginAbund(abund,el):
        if el not in Decay.LAs:
            print("Element "+el+" not in radiosotope list.")
        LAel=Decay.LAs[el]
        HLs=Decay.SSAGE/LAel
        return abund*(2**HLs)
        

    def calcHeatProd(el,conc, cp, dt):
        return (Decay.H0s[el]*conc/cp*dt)
        
    def calcDecayFrac(el,C0,time):
        return C0*np.exp(-Decay.LN2*time/Decay.LAs[el])
