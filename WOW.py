import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from decay import *
from thermalConduction import *
from aqStep import *
from matTrans import *
from rockEquil import *
from lookupProps import *

class cellTypes:
    UNDIFF=0
    MIXED=1
    FLUID=2
    ROCK=3
    ICE=4
    IGNORE=5

class gridCell:
    def __init__(self):
        self.Celltype = 0
        self.Time=0
        self.Temp=0
        self.Press=0
        self.Celltype = 0
        self.AqComp = dict()
        self.AqSpec = dict()
        self.RIComp = dict()
        self.IceComp= dict()
        self.RockComp=dict()
        self.RockPhases=[]
        self.RockPhaseDat=[]
        self.Top=0
        self.Bot=0
        self.Vol=0
        self.Mass=0
        self.RockMass=0
        self.IceMass=0
        self.AqMass=0
        self.Dens=0
        self.TCond=0
        self.Cp=0
        self.Porosity=0
        self.TempStep=50
        self.lastEquil=0
        self.lastEnthalpy=0
        self.latentHeat=0
        self.HtoAdd=0
        self.doUpdate=0

        
    def cell_decay(self,dtime):
        for i in self.RIComp:
            self.RIComp[i]=Decay.calcDecayFrac(i,self.RIComp[i],dtime)
        
    def cell_heat(self,dtime):
        DT = 0
        for i in self.RIComp:
            DT = DT + Decay.calcHeatProd(i, self.RIComp[i], self.Cp, dtime) #self.Vol, self.Cp, self.Dens)

        if self.latentHeat>0:
            self.latentHeat=self.latentHeat-(DT*self.Mass*self.Cp)
            if self.latentHeat<0:
                self.Temp = self.Temp-(self.latentHeat/self.Cp/self.Mass)
                self.latentHeat=0
        else:
            self.Temp+=DT
        
    def cell_equilibrate(self,k,perplex,eq36):
        self.RockMass=self.getRockMass()
        self.IceMass=self.getIceMass()
        self.AqMass=self.getWaterMass()
        if self.Celltype==cellTypes.FLUID:
            heat=callAqFreeze(self.AqComp, self.IceComp, self.Press, self.Temp)
            if eq36:
                #newComp,addH=aqStep.callAqEquil(self.AqComp, self.Press, self.Temp)
                newPrecip,newAq,addH=aqStep.callAqEquil({'CH4':0.0,'NH3':1.0,'CO':0.0,'CO2':0.0}, self.Press, self.Temp)
                self.AqSpec=newAq
                #heat+=addH
                self.doUpdate=0
            if heat>0:
                self.latentHeat=-heat
                self.Celltype=cellTypes.ICE
        if self.Celltype==cellTypes.ICE:
            #print('here')
            abc=1+1
            #callIceEquil(self.IceComp, self.AqComp, self.Press, self.Temp)
        if self.Celltype==cellTypes.ROCK and (self.Temp-self.lastEquil)>self.TempStep:
            if perplex:
                tempRC,tempRPD,tempRP,E,success=rockEquil(self.RockComp, self.Press, self.Temp)
                if success:
                    self.RockComp=tempRC
                    self.RockPhaseDat=tempRPD
                    self.RockPhases=tempRP
                    for key in self.RockComp:    
                        self.RockComp[key] *= self.RockMass #self.Mass-(self.getIceMass()+self.getWaterMass())
                    if not self.lastEnthalpy==0:
                        self.latentHeat+=(E-self.lastEnthalpy)*self.RockMass
                        if self.latentHeat<0:
                            self.HtoAdd=-self.latentHeat
                            self.latentHeat=0
                        self.lastEnthalpy=E
                    self.lastEquil=self.Temp
                    self.doUpdate=1
                #print("perplexed")
        if self.Celltype==cellTypes.UNDIFF:
            heat=0
            heat=callIceEquil(self.IceComp, self.AqComp, self.Press, self.Temp)
            if heat>0:
                self.Celltype=cellTypes.MIXED
                self.latentHeat=heat
                self.lastEquil=0
        if self.Celltype==cellTypes.MIXED:
            self.RockComp, self.AqComp, heat = aqStep.callAqRockEquil(self.RockComp, self.AqComp, self.Press, self.Temp,10)
        #if self.Celltype==cellTypes.IGNORE:
        #    self.AqComp={"H2O":1}
        return self.doUpdate
 
    def cell_updateProp(self):
        #RR=self.Temp-self.lastEquil #Is this robust? If I add a step between equil and update, probably not.
        if True: #self.doUpdate:
            self.TCond=LookupProps.calcThermalCond(self.Press,self.Temp,self.RockPhases,self.RockPhaseDat,self.IceComp,self.AqComp,self.Mass,self.RockComp)
            self.Cp = LookupProps.calcHeatCap(self.Press,self.Temp,self.RockPhases,self.RockPhaseDat,self.IceComp,self.AqComp,self.Mass,self.RockComp)
            self.Dens =LookupProps.calcDens(self.Press,self.Temp,self.RockPhases,self.RockPhaseDat,self.IceComp,self.AqComp,self.Mass,self.RockComp)
        if self.HtoAdd>0:
            self.Temp+=self.HtoAdd/self.Cp/self.Mass
            self.HtoAdd=0
        return

    def copy_cell(self,OldCell,time):
        self.Time=time
        self.Temp=OldCell.Temp
        self.Press=OldCell.Press
        self.Celltype = OldCell.Celltype
        self.AqComp = OldCell.AqComp.copy()
        self.AqSpec = OldCell.AqSpec.copy()
        self.RIComp = OldCell.RIComp.copy()
        self.IceComp= OldCell.IceComp.copy()
        self.RockComp=OldCell.RockComp.copy()
        self.RockPhases=OldCell.RockPhases.copy()
        self.RockPhaseDat=OldCell.RockPhaseDat.copy()
        self.Top=OldCell.Top
        self.Bot=OldCell.Bot
        self.Vol=OldCell.Vol
        self.Mass=OldCell.Mass
        self.AqMass=OldCell.AqMass
        self.IceMass=OldCell.IceMass
        self.RockMass=OldCell.RockMass
        self.Dens=OldCell.Dens
        self.TCond=OldCell.TCond
        self.Cp=OldCell.Cp
        self.Porosity=OldCell.Porosity
        self.TempStep=OldCell.TempStep
        self.lastEquil=OldCell.lastEquil
        self.latentHeat=OldCell.latentHeat
        
    def calc_vol(self):
        V1=4/3*np.pi*self.Top**3
        V2=4/3*np.pi*self.Bot**3
        self.Vol=V1-V2
    
    def calc_thickness(self,Bot):
        self.Bot=Bot 
        V2=4/3*np.pi*self.Bot**3
        self.Vol=self.Mass/self.Dens
        V1=self.Vol+V2
        T3 = V1/(4/3*np.pi)
        self.Top=T3**(1/3)
        return self.Top

    def calc_den_Mass(self):
        self.Mass=self.Vol*self.Dens

    def calc_Mass(self):
        self.Mass=self.getRockMass()+self.getIceMass()+self.getWaterMass()

    def getRockMass(self):
        sumR=0
        for i in self.RockComp:
            sumR=sumR+self.RockComp[i]
        return sumR
    
    def getIceMass(self):
        sumI=0
        for i in self.IceComp:
            sumI=sumI+self.IceComp[i]
        return sumI

    def getWaterMass(self):
        sumW=0
        for i in self.AqComp:
            sumW=sumW+self.AqComp[i]
        return sumW

    def removeWaterMass(self):
        WaterDict=self.AqComp.copy()
        self.AqComp={}
        return WaterDict
    
    def removeRockMass(self):
        rockDict=self.RockComp.copy()
        self.RockComp={}
        return rockDict

    def replaceMass(self,wcomp,rcomp):
        self.RockComp=rcomp
        self.AqComp=wcomp
        self.reclassify()

    def reclassify(self):
        if self.RockMass==0 and self.IceMass==0: #not self.RockComp:
            self.Celltype=cellTypes.FLUID
        elif self.IceMass==0 and self.AqMass==0: #not self.AqComp:
            self.Celltype=cellTypes.ROCK
        elif self.RockMass==0 and self.AqMass==0: #not self.AqComp:
            self.Celltype=cellTypes.ICE
        elif not self.Celltype==cellTypes.ICE:
            self.Celltype=cellTypes.IGNORE
            

class Planet:    
    def __init__(self,nr,startTime,endTime,Radius,maxnt,Tstep,Pstep):
        self.nr = nr
        self.nt = maxnt
        self.times=np.array([startTime])
        self.endTime=endTime
        self.boundTemp=0
        self.radii=np.linspace(0,Radius,nr+1) # Delete and replace with correct radii from each cell
        self.grid = [gridCell() for x in range(self.nr*self.nt)]
        self.grid = np.reshape(self.grid,[self.nt,self.nr])
        self.Radius = Radius
        self.PressStep=Pstep   
        self.Mass=0
        self.TempStep=Tstep
        self.eq36=0

    def timeStep(self,k):
        Neq36=0
        Yeq36=self.eq36
        self.eq36=0
        if k%100==0:
            print("Step "+str(k)+" out of "+str(self.nt))
            #print(self.Mass)
        for i in range(0,self.nr):
            if k>0:
                if i<self.nr-1:
                    self.grid[k,i].cell_heat(self.times[k]-self.times[k-1])
                self.grid[k,i].cell_decay(self.times[k]-self.times[k-1])
            if k==0:
                self.grid[k,i].cell_decay(self.times[k])
            if i%self.PressStep == 0:
                Neq36=self.grid[k,i].cell_equilibrate(k,1,Yeq36)
            else:
                Neq36=self.grid[k,i].cell_equilibrate(k,0,Yeq36)
            if Neq36:
                self.eq36=Neq36
        
        self.copyRockComp(k)
        
        for i in range(0,self.nr):    
            matTrans.extractPoreFluid(self.grid[k,i])
            self.grid[k,i].cell_updateProp()
        
        self.transportMaterial(k)
        
        for i in range(0,self.nr):    
            self.grid[k,i].cell_updateProp()
        
        self.transferHeat(k)
        self.updateBulk(k)
        
        for i in range(0,self.nr):
            self.grid[k+1,i].copy_cell(self.grid[k,i],self.times[k+1])
    
    def transferHeat(self,k):
        temp = np.array([self.grid[k,i].Temp for i in range(0,self.nr)])
        Cp =np.array([self.grid[k,i].Cp for i in range(0,self.nr)])
        Tcond = np.array([self.grid[k,i].TCond for i in range(0,self.nr)])
        rho = np.array([self.grid[k,i].Dens for i in range(0,self.nr)])
        Rs = np.array([self.grid[k,i].Bot for i in range(0,self.nr)]) # Is this right??? Use average of top and bottom?
        temp2, dtime = ThermalConduction.thermalCondStep(temp,Cp,Tcond,rho,Rs,k)
        self.times=np.append(self.times,self.times[k]+dtime)
        for i in range(0,self.nr):
            dT=temp2[i]-self.grid[k,i].Temp
            latentT=(self.grid[k,i].latentHeat/self.grid[k,i].Mass/self.grid[k,i].Cp)
            # this logic can be simplified:
            if dT>0:
                if latentT>dT:
                    self.grid[k,i].latentHeat=self.grid[k,i].latentHeat-(dT*self.grid[k,i].Mass*self.grid[k,i].Cp)
                    dT=0
                else:
                    dT=dT-latentT
                    self.grid[k,i].latentHeat=0
            if dT<0 and latentT<0:
                if latentT<dT:
                    self.grid[k,i].latentHeat=self.grid[k,i].latentHeat-(dT*self.grid[k,i].Mass*self.grid[k,i].Cp)
                    dT=0
                else:
                    dT=dT-latentT
                    self.grid[k,i].latentHeat=0
            self.grid[k,i].Temp=self.grid[k,i].Temp+dT
        
    def transportMaterial(self,time_step):
        topB=0
        bottomB=0
        iceBot=self.nr
        for i in range(0,self.nr):
            if self.grid[time_step,i].Celltype==cellTypes.MIXED:
                topB=i
            if self.grid[time_step,i].Celltype==cellTypes.ICE:
                break
                
            #if self.grid[time_step,i].Celltype==cellTypes.ROCK:
            #    bottomB=i
        if topB>bottomB: #and topB<IceBot:
            aqCopy, rockCopy, phaseDatCopy, phasesCopy=matTrans.differentiate(self.grid[time_step,:],topB,bottomB)
            for i in range(0,topB):
                #print(len(aqCopy))
                #print(i)
                self.grid[time_step,i].RockComp=rockCopy[i]
                self.grid[time_step,i].RockPhaseDat=phaseDatCopy[i]
                self.grid[time_step,i].RockPhases=phasesCopy[i]
                self.grid[time_step,i].AqComp=aqCopy[i]
            #print(bob)

    def copyRockComp(self,k):
        rockInd = 0
        swapp=0
        for i in range(self.nr-1,-1,-1):
            if self.grid[k,i].Celltype==cellTypes.ROCK:
                rockInd=i
                break
        for i in range(0,self.nr):
            if not len(self.grid[k,i].RockPhaseDat)==0:
                swapp=1
                break
        if swapp==1:
            lastStep=i
            for i in range(0,rockInd):
                if i%self.PressStep==0:
                    lastStep=i
                else:
                    self.grid[k,i].RockPhaseDat=self.grid[k,lastStep].RockPhaseDat.copy()
                    self.grid[k,i].RockPhases=self.grid[k,lastStep].RockPhases
                    sumC=0
                    for key in self.grid[k,lastStep].RockComp:
                        sumC += self.grid[k,lastStep].RockComp[key]
                    for key in self.grid[k,lastStep].RockComp:
                        self.grid[k,i].RockComp[key]=self.grid[k,lastStep].RockComp[key]/sumC*self.grid[k,i].RockMass
        
    def updateBulk(self,time_step):
        bot=0
        for i in range(0,self.nr):
            bot=self.grid[time_step,i].calc_thickness(bot)
            self.grid[time_step,i].calc_vol()
            self.grid[time_step,i].calc_Mass()
        self.calc_press(time_step)
       
 
    def initialize_comp(self,initIceComp,initRockComp, initRIComp, initTemp, initRho, initK, initCp):
        self.grid[0,0].IceComp=initIceComp
        self.grid[0,0].RockComp=initRockComp
        self.grid[0,0].RIComp=initRIComp
        self.grid[0,0].Temp=initTemp
        self.grid[0,0].Dens=initRho
        self.grid[0,0].TCond=initK
        self.grid[0,0].Cp=initCp
        self.grid[0,0].TempStep=self.TempStep
        self.boundTemp=initTemp
        for i in range(1,self.nr):
            self.grid[0,i].copy_cell(self.grid[0,0],self.times[0])
            self.grid[0,i].Bot=self.radii[i]
            self.grid[0,i].Top=self.radii[i+1]
            self.grid[0,i].calc_vol()
            self.grid[0,i].calc_den_Mass() ### FIXXXXXXXXXXXXX
            for j in self.grid[0,i].IceComp:
                self.grid[0,i].IceComp[j]=self.grid[0,i].IceComp[j]*self.grid[0,i].Mass
            for j in self.grid[0,i].RockComp:
                self.grid[0,i].RockComp[j]=self.grid[0,i].RockComp[j]*self.grid[0,i].Mass
        

        self.grid[0,0].Bot=0
        self.grid[0,0].Top=self.radii[1]
        self.grid[0,0].calc_vol()
        self.grid[0,0].calc_den_Mass()
        for j in self.grid[0,0].IceComp:
            self.grid[0,0].IceComp[j]=self.grid[0,0].IceComp[j]*self.grid[0,0].Mass
        for j in self.grid[0,0].RockComp:
            self.grid[0,0].RockComp[j]=self.grid[0,0].RockComp[j]*self.grid[0,0].Mass
        self.Mass=self.calc_mass(0)
        self.calc_press(0)

 
    def plotAttribute(self,att,tit):
        matplotlib.rcParams.update({'font.size': 22})
        data = np.zeros((self.nt,self.nr))
        for i in range(self.nr):
            for j in range(self.nt):
                data[j][i]=getattr(self.grid[j][i], att)
        plt.figure(figsize=(10,6))
        #print(np.shape(self.radii[0:-1]))
        #print(np.shape(self.times))
        #print(np.shape(self.grid))
        plt.contourf(self.times/Decay.YR/10**6,self.radii[0:-1]/1000,np.transpose(data),100) # Fix to use correct radii
        plt.colorbar()
        plt.contour(self.times/Decay.YR/10**6,self.radii[0:-1]/1000,np.transpose(data),15, colors='k')
        plt.xlabel('Time (Myr)')
        plt.ylabel('Radius (km)')
        plt.title(tit)
        
    def plotAttributeLine(self,att,tit,radius,unit=1):
        matplotlib.rcParams.update({'font.size': 22})
        data = np.zeros(self.nt)
        for j in range(self.nt):
            data[j]=getattr(self.grid[j][radius], att)
        plt.figure(figsize=(10,6))
        plt.plot(self.times/Decay.YR/10**6,data*unit,100) # Fix to use correct radii
        plt.xlabel('Time (Myr)')
        plt.ylabel(tit)
        #print(data)
    
    def plotDictAttribute(self,att,key,tit,unit=1):
        matplotlib.rcParams.update({'font.size': 22})
        data = np.zeros((self.nt,self.nr))
        for i in range(self.nr):
            for j in range(self.nt):
                tempItem=getattr(self.grid[j][i], att)
                if key in tempItem:
                    data[j][i]=tempItem[key]
        plt.figure(figsize=(10,6))
        plt.contourf(self.times/Decay.YR/10**6,self.radii[0:-1]/1000,np.transpose(data)*unit,100) # Fix to use correct radii
        plt.colorbar()
        plt.contour(self.times/Decay.YR/10**6,self.radii[0:-1]/1000,np.transpose(data)*unit,15, colors='k')
        plt.xlabel('Time (Myr)')
        plt.ylabel('Radius (km)')
        plt.title(tit)
        
    def plotDictAttributeMassScaled(self,att,key,tit,vmin,vmax,unit=1):
        matplotlib.rcParams.update({'font.size': 22})
        data = np.zeros((self.nt,self.nr))
        for i in range(self.nr):
            for j in range(self.nt):
                tempItem=getattr(self.grid[j][i], att)
                if key in tempItem:
                    data[j][i]=tempItem[key]/self.grid[j][i].Mass
        plt.figure(figsize=(10,6))
        plt.contourf(self.times/Decay.YR/10**6,self.radii[0:-1]/1000,np.transpose(data)*unit,100, vmin=vmin,vmax=vmax) # Fix to use correct radii
        plt.colorbar()
        plt.contour(self.times/Decay.YR/10**6,self.radii[0:-1]/1000,np.transpose(data)*unit,15, colors='k')
        plt.xlabel('Time (Myr)')
        plt.ylabel('Radius (km)')
        plt.title(tit)
    
    def plotDictAttributeLine(self,att,key,tit,radius,unit=1):
        matplotlib.rcParams.update({'font.size': 22})
        data = np.zeros(self.nt)
        for j in range(self.nt):
            tempItem=getattr(self.grid[j][radius], att)
            data[j]=0
            if key in tempItem:
                data[j]=tempItem[key]
        plt.figure(figsize=(10,6))
        plt.plot(self.times/Decay.YR/10**6,data*unit,100) # Fix to use correct radii
        plt.xlabel('Time (Myr)')
        plt.ylabel(tit)
        #print(data)
        return plt.gca()
   
    def plotPhaseAssemblage(self,radius):
        pP=dict()
        for j in range(self.nt):
            phases=self.grid[j][radius].RockPhases
            for i in range(0,len(phases)):
                if not phases[i]=="Bulk_rs":
                    if phases[i] in pP:
                        pP[phases[i]][j]=self.grid[j][radius].RockPhaseDat[i]["wt%"]
                    else:
                        pP[phases[i]]=np.zeros(self.nt)
                        pP[phases[i]][j]=self.grid[j][radius].RockPhaseDat[i]["wt%"]
        
        #x_values = list(pP.keys())
        #y_values_list = list(pP.values())
   
        #print(pP.values())
        #print(pP.keys())

        plt.figure(figsize=(10,6))
        plt.stackplot(self.times/Decay.YR/10**6,pP.values(),labels=pP.keys()) 
        #plt.stackplot(x_values, y_values_list, labels=list(pP.keys()))
        plt.legend(loc='upper left')
        plt.xlabel('Time (Myr)')
        plt.ylabel('Wt %')
        plt.title('Phase assemblage at %0.2f kms deep'%(self.radii[-1]-self.radii[radius]))
        plt.show()


    def runModel(self):
        i=0
        while self.times[-1] < self.endTime and i<self.nt-1: #i in range(0,nt-1):
            #print(self.times[i])
            self.timeStep(i)
            i=i+1

    def calc_press(self,time_index):
        self.Mass=self.calc_mass(time_index)
        M=self.Mass
        dP=0
        for i in range(self.nr-1,0,-1):
            G = 6.67*10**-11
            R = self.grid[time_index][i].Bot
            M = M-self.grid[time_index][i].Mass
            g = G*M/R**2
            dh = self.grid[time_index][i].Top-R
            dens=self.grid[time_index][i].Dens
            dP = dP+dens*g*dh
            self.grid[time_index][i].Press=dP
            #print(R)
            #print(M)
            #print(g)
            #print(dh)
            #print(dens)
        self.grid[time_index][0].Press=dP

    def calc_mass(self,time_index):
        sumM=0
        for i in range(0,self.nr):
            sumM = sumM+self.grid[time_index][i].Mass
        return sumM



######### User helper functions


def calcOriginAbundance(abund,el):
    return Decay.calcOriginAbund(abund,el)



######### To be omitted


def callAqFreeze(waterDict, iceDict, press, temp):
    heat = 0
    IceMeltHeat=334*1000 #J/kg
    if temp<273.15:
        if 'H2O' in waterDict:
            iceDict['H2O']=waterDict['H2O']
            heat=waterDict["H2O"]*IceMeltHeat
            waterDict['H2O']=0
    return heat
    
def callIceEquil(IceDict, waterDict, press, temp):
    heat = 0
    IceMeltHeat=334*1000 #J/kg
    if temp>273.15:
        waterDict["H2O"]=IceDict["H2O"]
        IceDict["H2O"]=0
        heat=waterDict["H2O"]*IceMeltHeat
    return heat
        














