import numpy as np

class ThermalConduction:
    def calcStep(Cp, K, rho, DR, timeC): # WRONG!
        if timeC<5*3.1558e+13: # 3myr
            return 3.1558e+12 # debugging test
        alpha = K/(Cp*rho)
        DT = np.min(DR**2/alpha/3)
        return DT

    def thermalCondStep(T,Cp,K,rho,Rads,k,timeC):
        nc = len(T)

        C = Rads
        DR = np.diff(C,append=0)
        DR[-1]=DR[-2]

        DT = ThermalConduction.calcStep(Cp,K,rho,DR,timeC)

        D=np.zeros((nc,3))

        for i in range(0,nc):
            for j in range(0,3):
                if i==0:
                    D[0,j]=K[0]/(rho[i]*Cp[0])
                elif i==nc-1:
                    D[nc-1,j]=K[nc-1]/(rho[i]*Cp[nc-1])
                else:
                    D[i,j]=(K[i]+K[i+(j-1)])/(rho[i]*Cp[i])/2.0

        dn = D[:,1]
        dp = D[:,2]
        dm = D[:,0]

        diagU = np.zeros(nc)
        diagC = np.zeros(nc)
        diagL = np.zeros(nc)

        for i in range(1,nc-1):
            diagL[i] = DT/(2.0*DR[i]*DR[i+1])*(dn[i]/(i)-dm[i])
            diagC[i] = 0.5+DT/(2.0*DR[i]*DR[i+1])*(dp[i]+dm[i])
            diagU[i] = -DT/(2.0*DR[i]*DR[i+1])*(dn[i]/(i)+dp[i])

        diagC[0] = 0.5 + DT/(DR[0]**2)*dn[0]
        diagU[0] = -DT/(DR[0]**2)*dn[0]

        # Calc source
        Hadd = np.zeros(nc)
        for i in range(0,nc-1):
            Hadd[i] = 0.5*T[i] #+DT*H[i]/(2*Cp[i])
        Hadd[nc-2]=Hadd[nc-2]+DT/(2*DR[nc-2]**2)*(dn[nc-2]/(nc-2)+dp[nc-2])*T[nc-1]

        # Solve
        U = np.zeros(nc)
        GAM = np.zeros(nc)

        BET=diagC[0]
        U[0] = Hadd[0]/BET
        U[-1]=T[-1]
        #print(k)
        #print(T[-1])
        #print(T[-2])

        for i in range(1,nc-1):
            GAM[i]=diagU[i-1]/BET
            BET=diagC[i]-diagL[i]*GAM[i]
            U[i] = (Hadd[i]-diagL[i]*U[i-1])/BET
        for i in range(nc-2,-1,-1):
            U[i]=U[i]-GAM[i+1]*U[i+1]

        return U,DT

