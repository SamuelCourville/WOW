import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from molmass import Formula

from decay import *
from thermalConduction import *
from aqStep import *
from matTrans import *
from rockEquil import *
from lookupProps import *
from plotHelper import *

# System to define cell types
class cellTypes:
    UNDIFF=0  # Use for cells that are primordial frozen accreted material
    MIXED=1   # Use for cells that are at compositional boundaries (e.g., rocky-mantle/ocean interface or ice/ocean interface)
    FLUID=2   # Use for cells that are fluid (e.g., ocean)
    ROCK=3    # Use for cells that are part of the rocky interior (may have fluid in pore space in future implementations).
    ICE=4     # Use for cells that are part of the outer ice shell
    IGNORE=5  # Use for testing or error catching





##############################################
#
#  Grid cell class
#
##############################################

# Class for grid cells
class gridCell:
    def __init__(self):
        '''
        A grid cell is one radius increment of the full planet/moon/asteroid model. A full model is composed of many grid cells.
        A grid cell retains the same mass throughout the entire simulation, though its volume may change as its composition and density does.
        The grid cell contains all cell specific properties including composition and physical properties.
        '''

        # Basic cell definition and states
        self.Celltype = 0     # See cellType class above
        self.ind = 0          # what position the cell is in in the grid
        self.Time=0           # The current time step in the model (is this redundant? Not every cell needs its own time counter)
        self.Top=0            # The top bound of the cell, in meters, measured from the core.
        self.Bot=0            # The bottom bound of the cell, in meters, measured from the core.
        self.Vol=0            # The volume of the cell
        self.Mass=0           # The mass of the cell (This is defined to stay constant once initialized)
        self.Temp=0           # The current temperature of the cell
        self.MaxOrgTemp=0     # The maximum temperature thet the organics have reached in the cell NEEDS TO BE UPGRADED! BAD IDEA
        self.Press=0          # The current pressure within the cell
        # Cell properties
        self.Dens=0           # The current density of the cell
        self.DenCo=[0,0,0]    # ADD NOTE
        self.TCond=0          # The thermal conductivity of the cell
        self.TCondCo=[0,0,0]  # ADD NOTE
        self.Cp=0             # The heat capacity of the cell
        self.CpCo=[0,0,0]     # ADD NOTE
        self.Porosity=0       # The porosity of the cell (NOT IMPLEMENTED YET)

        # Cell aqueous composition
        self.pH=0             # pH of the fluid component  
        self.pH2=0            # pH2 of the fluid component  
        self.fO2=0            # fO2 (e.g., redox state) of the fluid component  
        self.AqComp = dict()  # A dictionary to store the elemental composition of the aqueous component
        self.AqSpec = dict()  # A dictionary to store the aqueous speciation of the aqueous composition
        self.AqAct = dict()  # A dictionary to store the aqueous speciation of the aqueous composition
        self.AqMin = dict()   # A dictionary to store the precipitated minerals of the aqueous composition
        self.AqGas = dict()   # A dictionary to store the aqueous gasses of the aqueous composition
        self.AqMass=0         # A helper variable to store just the mass of the aqueous component
        
        # Cell radiosotopes
        self.RIComp = dict()  # A dictionary to store radioisotope abundances

        # Cell ice composition
        self.IceComp= dict()  # A dictionary to store the composition of the ice component
        self.IceMass=0        # A helper variable to store the mass of the ice component
        
        # Cell rock compositions
        self.RockComp=dict()  # A dictionary to store elemental abundances of the rock composition
        self.RockPhases=[]    # An array to store the list of mineral phases
        self.RockPhaseDat=[]  # An array to store dictionaries of properties for each mineral phase
        self.RockMass=0       # A helper variable to store the mass of just the rock component 
        #self.OrigOrgMass=0    # what wt% of the cell was originally IOM
        self.OrgComp=dict()   # A dictionary to store elemental abundances of the organics in the rock composition
        
        # Flags and condition variables
        self.TempStep=50      # How many degrees to change before mineralogy will be equilibrated
        self.lastEquil=0      # Last quilibration temperature
        self.lastEnthalpy=0   # Last enthalpy of formation for the mineral phase assemblage (J/kg)
        self.latentHeat=0     # A buffer to store latent heat that must be resolved (J)
        self.latentFreeze=0   # A buffer to store latent heat from freezing (J)
        self.HtoAdd=0         # Extra heat to add (K)
        self.doUpdate=0       # Flag to decide whether or not to reequilibrate aqueous component

        # Currently, when thermodynamic equilibrium predicts a change in composition, it saves the ammount of heat required, or produced,
        # into the latent heat variable. Temperature of the cell cannot be increased until a latent heat deficit is accounted for. 
        # This method is backwards. The latent heat should be addressed in concurance with equilibrium, not afterwards. But that's tough to do.
        
    def cell_decay(self,dtime):
        '''
           Updates the abundances of radioactive isotopes after decaying for a set interval of time, dtime. 
        '''
        for i in self.RIComp:
            self.RIComp[i]=Decay.calcDecayFrac(i,self.RIComp[i],dtime)
        
    def cell_heat(self,dtime):
        '''
           Calculates the heat produced from decaying radioactive isotopes over the time interval, dtime. Updates the temperature of the cell.
        '''
        DT = 0   # change in temperature
        for i in self.RIComp:
            DT = DT + Decay.calcHeatProd(i, self.RIComp[i], self.Cp, dtime)

        # applies heat to the latent heat deficit buffer first, then changes the cell temperature if not zero.
        if self.latentHeat>0:
            self.latentHeat=self.latentHeat-(DT*self.Mass*self.Cp)
            if self.latentHeat<0:
                self.Temp = self.Temp-(self.latentHeat/self.Cp/self.Mass)
                self.latentHeat=0
        else:
            self.Temp+=DT
        
    def cell_equilibrate(self,k,perplex,eq36,name):
        '''
           Equilibrates the composition of the cell depending on its temperature and pressure.
           k: the current time step. Used for debugging I think? Not needed anymore.
           perplex: A boolean flag to decide whether to reequilibrate the mineral assemblage with perplex on this iteration or not.
           eq36: A boolean flag to decide whether to reequilibrate the fluid with eq36 on this iteration.

           returns doUpate: a flag to decide whether eq36 needs to be run on the next iteration or not (decided by whether perplex ran)
        '''
        # Gets cell type
        ct = self.Celltype
        
        # Updates rock, ice, and fluid masses before equilibration steps
        self.RockMass=self.getRockMass()  
        self.IceMass=self.getIceMass()    
        self.AqMass=self.getWaterMass()   

        iceImp={}
        # If statements to decide which type of equilibration is needed based on what type of composition the cell has.

        # Fluid cell uses eq36 equilibration (NOT COMPLETE)
        if ct==cellTypes.FLUID:
            self.AqComp = aqStep.refactor_water(self.AqComp) # removes O and H and changes it to H2O in dictionary
            heat=0
            heat2,frozen,newIceComp,iceImp=aqStep.freeze(self.AqComp, self.Press, self.Temp, name,self.ind)
            if frozen:
                self.IceComp=newIceComp
                self.AqComp={}
                self.AqSpec={}
                self.AqAct={}
                self.pH=np.nan
                self.pH2=np.nan
                self.AqGas={}
                self.AqMin={}
                self.Celltype=cellTypes.ICE
            if "H2O" in self.AqComp and not frozen and eq36: # I changed this to exclude temperature.
                newpH, newpH2, newMin,newAq, newAct, newGas,addH=aqStep.callAqEquil(self.AqComp, self.Press, self.Temp,name,self.ind)
                self.AqSpec=newAq
                self.AqAct = newAct
                self.pH=newpH
                self.pH2=newpH2
                self.AqGas=newGas
                self.AqMin=newMin
                heat+=addH
                self.doUpdate=0
            #heat+=heat2
            self.latentFreeze=heat2
            #if heat>0:
            #    self.latentHeat=-heat
 
        # Ice shell uses FREZCHEM equilibration (NOT YET IMPLEMENTED!!!!!!)
        if ct==cellTypes.ICE:
            heat=callIceEquil(self.IceComp, self.AqComp, self.Press, self.Temp)
            if heat>0:
                self.Celltype=cellTypes.FLUID
                self.latentHeat=heat
                self.lastEquil=0

        # Rock cell uses perplex equilibration
        if ct==cellTypes.ROCK and (self.Temp-self.lastEquil)>self.TempStep:
            if perplex:
                tempRC,tempRPD,tempRP,E,NOT,vols,success=rockEquil(self.RockComp, self.Press, self.Temp, self.MaxOrgTemp)
                if success:
                    self.MaxOrgTemp=NOT
                    self.RockComp=tempRC
                    self.RockPhaseDat=tempRPD
                    self.RockPhases=tempRP
                    self.OrgComp=vols
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
                else:
                    self.lastEquil=self.Temp # DEBUGGING HACK. GET RID OF
                #print("perplexed")
        
        # Undifferentiated cell does not use thermodynamics to equilibrate. Melts at 273K. This should be improved.
        if ct==cellTypes.UNDIFF:
            heat=0
            heat=callIceEquil(self.IceComp, self.AqComp, self.Press, self.Temp)
            if heat>0:
                self.Celltype=cellTypes.MIXED
                self.latentHeat=heat
                self.lastEquil=0

        # Mixed cell to handle interface condtions (NOT COMPLETE) 
        if ct==cellTypes.MIXED:
            self.RockComp, self.AqComp, heat = aqStep.callAqRockEquil(self.RockComp, self.AqComp, self.Press, self.Temp,10)
        #if self.Celltype==cellTypes.IGNORE:
        #    self.AqComp={"H2O":1}

        # Return necessary flags for future equilibration steps
        return self.doUpdate, iceImp
 
    def cell_updateProp(self, ks, rhos, cps, crossRef):
        '''
           Updates the physical and thermal properties of the cell, to be done after an equilibration step.
        '''
        if True: #self.doUpdate:
            self.TCond=LookupProps.calcThermalCond(ks,crossRef,self.Press,self.Temp,self.RockPhases,self.RockPhaseDat,self.IceComp,self.AqComp,self.Mass,self.RockComp)
            self.Cp= LookupProps.calcHeatCap(cps,crossRef,self.Press,self.Temp,self.RockPhases,self.RockPhaseDat,self.IceComp,self.AqComp,self.Mass,self.RockComp)
            self.Dens=LookupProps.calcDens(rhos,crossRef,self.Press,self.Temp,self.RockPhases,self.RockPhaseDat,self.IceComp,self.AqComp,self.Mass,self.RockComp)
        else: # attempt to save time, but no way for it to work
            self.TCond=LookupProps.TCondFun(self.TCondCo,self.Press,self.Temp)
            self.Cp=LookupProps.HeatCapFunc(self.CpCo,self.Press,self.Temp)
            self.DenCo=LookupProps.rhoFunc(self.DenCo,self.Press,self.Temp)
        if self.HtoAdd>0:
            self.Temp+=self.HtoAdd/self.Cp/self.Mass
            self.HtoAdd=0
        return

    def copy_cell(self,OldCell,time):
        '''
           Used to copy the compoistion and properties of one cell to a new cell. 
           Inputs: 
              time: the current time step
              OldCell: the previous cell to be copied
        '''
        self.Time=time
        self.ind=OldCell.ind
        self.Temp=OldCell.Temp
        self.MaxOrgTemp=OldCell.MaxOrgTemp
        self.Press=OldCell.Press
        self.Celltype = OldCell.Celltype
        self.AqComp = OldCell.AqComp.copy()
        self.AqSpec = OldCell.AqSpec.copy()
        self.AqAct = OldCell.AqAct.copy()
        self.AqMin = OldCell.AqMin.copy()
        self.AqGas = OldCell.AqGas.copy()
        self.pH=OldCell.pH
        self.pH2=OldCell.pH2
        self.fO2=OldCell.fO2
        self.RIComp = OldCell.RIComp.copy()
        self.IceComp= OldCell.IceComp.copy()
        self.RockComp=OldCell.RockComp.copy()
        #self.OrgComp=OldCell.OrgComp.copy()
        #self.OrigOrgMass=OldCell.OrigOrgMass
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
        self.DenCo=OldCell.DenCo
        self.TCond=OldCell.TCond
        self.TCondCo=OldCell.TCondCo
        self.Cp=OldCell.Cp
        self.CpCo=OldCell.CpCo
        self.Porosity=OldCell.Porosity
        self.TempStep=OldCell.TempStep
        self.lastEquil=OldCell.lastEquil
        self.latentHeat=OldCell.latentHeat
        self.latentFreeze = OldCell.latentFreeze
        
    def calc_vol(self):
        '''
           Calculates the volume of the cell from the upper and lower radii of the cell. 
        '''
        V1=4/3*np.pi*self.Top**3
        V2=4/3*np.pi*self.Bot**3
        self.Vol=V1-V2
    
    def calc_thickness(self,Bot):
        '''
           Calculates the thckness of the cell from the mass, density, and radial extent of the bottom of the cell. 
           Designed to be done in a loop over all grid cells, solving from the core to surface, passing the top radius of one cell
           to be the bottom radius of the next cell.
           
           Input: Bot, the bottom of the cell, measured as distance from the core in meters.
           Output: Top, the top of the cell in meters.
        '''
        self.Bot=Bot 
        V2=4/3*np.pi*self.Bot**3
        self.Vol=self.Mass/self.Dens
        V1=self.Vol+V2
        T3 = V1/(4/3*np.pi)
        self.Top=T3**(1/3)
        return self.Top

    # Is this funciton still necessary?
    def calc_den_Mass(self):
        '''
           Calculates the mass of the cell from the cells volume and density.
        '''
        self.Mass=self.Vol*self.Dens

    # Is this function necesssary or a good idea? The mass of each cell should never change.
    def calc_Mass(self):
        '''
           Calculates the mass of the cell from the sum of the masses of the rock, ice, and fluid components
        '''
        self.Mass=self.getRockMass()+self.getIceMass()+self.getWaterMass()

    def getRockMass(self):
        '''
           Adds up the masses of all the mineral components
        '''
        sumR=0
        for i in self.RockComp:
            sumR=sumR+self.RockComp[i]
        return sumR
    
    def getIceMass(self):
        '''
           Adds up the masses of all the ice components
        '''
        sumI=0
        for i in self.IceComp:
            sumI=sumI+self.IceComp[i]
        return sumI

    def getWaterMass(self):
        '''
           Adds up the masses of all the water components
        '''
        sumW=0
        for i in self.AqComp:
            sumW=sumW+self.AqComp[i]
        return sumW

    def removeWaterMass(self):
        '''
           Sets aqueous composition to empty, and returns existing composition
        '''
        WaterDict=self.AqComp.copy()
        self.AqComp={}
        return WaterDict
    
    def removeRockMass(self):
        '''
           Sets rock composition to empty, and returns existing composition
        '''
        rockDict=self.RockComp.copy()
        self.RockComp={}
        return rockDict

    # Is this still necessary
    def replaceMass(self,wcomp,rcomp):
        '''
           update the rock and water dictionaries with new compositions
        '''
        self.RockComp=rcomp
        self.AqComp=wcomp
        self.reclassify()
     
    # This funciton should be improved. If its even necessary?
    def reclassify(self):
        '''
           Checks the cell's compoisiton and updates the cell type accordingly. 
        '''
        if self.getRockMass()==0 and self.getIceMass()==0: #not self.RockComp:
            self.Celltype=cellTypes.FLUID
        elif self.getIceMass()==0 and self.getWaterMass()==0: #not self.AqComp:
            self.Celltype=cellTypes.ROCK
        elif self.getRockMass()==0 and self.getWaterMass()==0: #not self.AqComp:
            self.Celltype=cellTypes.ICE
        elif not self.getRockMass()==0 and not self.getWaterMass()==0:
            self.Celltype=cellTypes.MIXED
        elif not self.getRockMass()==0 and not self.getIceMass()==0:
            self.Celltype=cellTypes.UNDIFF
        elif not self.Celltype==cellTypes.ICE:
            self.Celltype=cellTypes.IGNORE
            

############################################################
#
#   Planet Class
#
############################################################

class Planet:    
    def __init__(self,nr,startTime,endTime,Radius,maxnt,Tstep,Pstep):
        self.name = ""          # what the model name is
        self.nr = nr            # number of grid cells in the model
        self.nt = maxnt         # Maximum number of time steps to take
        self.times=np.array([startTime]) # Array to keep track of time steps
        self.endTime=endTime    # endTime, if maxnt is not reached
        self.boundTemp=0        # surface temperature
        self.radii=np.linspace(0,Radius,nr+1) # Start each cell with the same thickness
        self.grid = [gridCell() for x in range(self.nr*self.nt)] # Initialize an array of grid cells with nr*nt cells
        self.grid = np.reshape(self.grid,[self.nt,self.nr]) # reshape array for more intuitive indexing
        self.Radius = Radius    # Radius of planet
        self.PressStep=Pstep    # Spacing to determine which cells equilibrate. Equilibrate every Pstep cells, and interpolate between for speed
        self.Mass=0             # Total mass of planet. Should stay the same
        self.TempStep=Tstep     # How many degrees K does a cell have to change before instructing the cell to reequilibrate.
        self.eq36=0             # Flag to decide whether or not to do an ocean equilibration step.
        self.extracts=[{}]*maxnt# Array of dictionaries to store extracted fluid at each time step
        self.kTable= LookupProps.loadk()
        self.rhoTable=LookupProps.loadrho()
        self.cpTable=LookupProps.loadcp()
        self.perpTable=LookupProps.loadperp()

    def timeStep(self,k):
        '''
           Backbone of modeling software. Calls funcitons to calculate thermal and chemical evolution of each cell.
        '''
        Neq36=0          # Flag to determine whether or not to call EQ36 for aqueous equilibration on this time step
        Yeq36=self.eq36  # Flag to determine whether or not to call EQ36 for aqueous equilibration on this time step 
        self.eq36=0      # Reset internal EQ36 flag
        
        if k%10==0:     # Print a statement every 100 time steps so the user can be sure the program is running. 
            print("Step "+str(k)+" out of "+str(self.nt))

        iceImp={}
        #### CALCUALTE DECAY HEAT
        # Loop through all cells
        for i in range(0,self.nr):
            # for first time step only. Calculate how much decay occurs before model starts. 
            if k==0:
                self.grid[k,i].cell_decay(self.times[k])
            # for all time steps past the first step, calculate decay and decay heat for each time step.
            if k>0:
                # do not add heat to cell on surface
                if i<self.nr-1:
                    self.grid[k,i].cell_heat(self.times[k]-self.times[k-1])
                self.grid[k,i].cell_decay(self.times[k]-self.times[k-1])

            #### EQUILIBRATE CELLS
            # If on one of the psteps, do a perplex equilibration calculation
            if i%self.PressStep == 0:
                Neq36,iceImpT=self.grid[k,i].cell_equilibrate(k,1,Yeq36,self.name)
            else:
                Neq36,iceImpT=self.grid[k,i].cell_equilibrate(k,0,Yeq36,self.name)
            # If a perplex equilibration step was just completed, set the next step to conduct an ocean eqilibration step.
            if Neq36:
                self.eq36=Neq36
            for ii in iceImpT:
                if ii in iceImp:
                    iceImp[ii]+=iceImpT[ii]
                else:
                    iceImp[ii]=iceImpT[ii]
            self.grid[k, i].reclassify()
        
        # For all cells that aren't on the perplex equilibration grid, copy values from the nearest equilibration step.
        self.copyRockComp(k)

        #### TRANSPORT MATERIAL
        self.extracts[k]=dict() 
        # Extract fluid generated by the Perplex equilibration, and then update the properties of the cell
        for i in range(0,self.nr):    
            updateNeeded, extracted=matTrans.extractPoreFluid(self.grid[k,i])  ###### Check that this works!!!!!!!!!!
            for i in extracted:
                if i in self.extracts:
                    self.extracts[k][i]+=extracted[i]
                else:
                    self.extracts[k][i]=extracted[i]
            #if updateNeeded:
            #    self.grid[k,i].cell_updateProp()  # Is this necessary? Why?
        # Differentiate the body if necessary (e.g., water from rock)
        self.transportMaterial(k,iceImp)

        #### UPDATE CELL PROPERTIES
        # Update cell properties after differentiation
        if True: #updateNeeded:
            for i in range(0,self.nr):    
                self.grid[k,i].cell_updateProp(self.kTable,self.rhoTable,self.cpTable,self.perpTable)

        #### TRANSFER HEAT
        # Do a thermal heat conduction step
        self.transferHeat(k)

        #### UPDATE BULK PROPERTIES
        # Update the bulk properties of the full grid
        self.updateBulk(k)
        
        # Copy cells into next time step
        for i in range(0,self.nr):
            self.grid[k+1,i].copy_cell(self.grid[k,i],self.times[k+1])
    
    def transferHeat(self,k):
        '''
           Does a thermal heat conduction step, dynamically calculates what the next time step should be. 
        '''
        temp = np.array([self.grid[k,i].Temp for i in range(0,self.nr)])   # Grabs an array of the current temps of all the cells
        Cp =np.array([self.grid[k,i].Cp for i in range(0,self.nr)])        # Grabs an array of the current Cp of all cells
        Tcond = np.array([self.grid[k,i].TCond for i in range(0,self.nr)]) # Grabs an array of the current thermal conductivity of all cells
        rho = np.array([self.grid[k,i].Dens for i in range(0,self.nr)])    # Grabs an array of the density of all cells
        Rs = np.array([self.grid[k,i].Bot for i in range(0,self.nr)])      # Is this right??? Use average of top and bottom?
        
        # Do a thermal conduction step and calculate the time step
        temp2, dtime = ThermalConduction.thermalCondStep(temp,Cp,Tcond,rho,Rs,k,self.times[k])
        
        # determine how to dstribute the temperature change from thermal conduction to account for latent heat.
        # This is not a good way to do it and should be fixed in the future. Currently, if heat is moved into a cell, the heat is first used 
        # to take away heat in the latent heat buffer.
        self.times=np.append(self.times,self.times[k]+dtime)
        for i in range(0,self.nr):
            dT=temp2[i]-self.grid[k,i].Temp
            latentT=(self.grid[k,i].latentHeat/self.grid[k,i].Mass/self.grid[k,i].Cp)
            latentF = (self.grid[k, i].latentFreeze / self.grid[k, i].Mass / self.grid[k, i].Cp)
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
            if dT<0 and latentF>0:
                if latentF>abs(dT):
                    self.grid[k, i].latentFreeze=self.grid[k, i].latentFreeze+(dT * self.grid[k, i].Mass * self.grid[k, i].Cp)
                    #print(dT)
                    dT=0
                    #print(dT)
                    #print(latentF)
                else:
                    #print(dT)
                    dT=dT+latentF
                    self.grid[k, i].latentFreeze = 0
                    #print(dT)
                    #print(latentF)
            self.grid[k,i].Temp=self.grid[k,i].Temp+dT
            if self.grid[k,i].Celltype==cellTypes.FLUID and self.grid[k,i].Temp>273.15:
                lostH=self.grid[k,i].Temp-273.15
                self.grid[k,i+1].Temp+=lostH # Send heat upwards
                self.grid[k, i].Temp=273.15
    
    
    def transportMaterial(self,time_step,iceImp):
        '''
           This function decides if material needs to move. E.g., if fluid is beneath rock, differentiate it.
        '''
        topB=self.nr            # Top of the rock and ocean
        bottomB=0         # Bottom of the rock
        #iceBot=self.nr
        #for i in range(0,self.nr):
        #    if self.grid[time_step,i].Celltype==cellTypes.UNDIFF:  # THIS MUST BE WRONG!!!!!!!!!!!
        #        topB=i    # find the top of the rocky mantle
        #    if self.grid[time_step,i].Celltype==cellTypes.ICE:
        #        break     # Stop looking if we reach ice
                
            # Not necessary anymore? What if we eventually get to rock/metal differentiation?
            #if self.grid[time_step,i].Celltype==cellTypes.ROCK:
            #    bottomB=i
        
        # differentiate all cells beneath the ocean. Does this always need to be done? Only should be done after an equilibrium change?
        if topB>bottomB: #and topB<IceBot:
            aqCopy, rockCopy, phaseDatCopy, phasesCopy=matTrans.differentiate(self.grid[time_step,:],topB,bottomB,iceImp)
            for i in range(0,topB): 
                self.grid[time_step,i].RockComp=rockCopy[i]
                self.grid[time_step,i].RockPhaseDat=phaseDatCopy[i]
                self.grid[time_step,i].RockPhases=phasesCopy[i]
                self.grid[time_step,i].AqComp=aqCopy[i]
                #self.grid[time_step,i].OrgComp=orgCopy[i]


    def copyRockComp(self,k):
        '''
           Interpolates rock composition from cells with Perplex equilibration to cells without.
        '''
        # Find index of highest rock cell
        rockInd = 0
        for i in range(self.nr-1,-1,-1):
            if self.grid[k,i].Celltype==cellTypes.ROCK:
                rockInd=i
                break
        # Find cells with data from Perplex
        swapp=0
        for i in range(0,self.nr):
            if not len(self.grid[k,i].RockPhaseDat)==0:
                swapp=1
                break
        # If there are cells with perplex data, continue.
        if swapp==1:
            lastStep=i
            for i in range(0,rockInd):   # Loop over only the rock cells
                if i%self.PressStep==0:  # Determine if the cell is on the perplex grid
                    lastStep=i
                else:                    # if not, copy data from the nearest cell below with perplex data
                    self.grid[k,i].RockPhaseDat=self.grid[k,lastStep].RockPhaseDat.copy()
                    self.grid[k,i].RockPhases=self.grid[k,lastStep].RockPhases
                    #self.grid[k,i].OrgComp=self.grid[k,lastStep].OrgComp.copy()
                    #self.grid[k,i].MaxOrgTemp=self.grid[k,lastStep].MaxOrgTemp
                    self.grid[k,i].lastEquil=self.grid[k,lastStep].lastEquil
                    
                    # Rescale the copied rock comosition to match the mass of the cell.
                    sumC=0
                    for key in self.grid[k,lastStep].RockComp:
                        sumC += self.grid[k,lastStep].RockComp[key]
                    #print(sumC)
                    #print(lastStep)
                    #print(self.grid[k,lastStep].RockComp)
                    #print(self.grid[k, lastStep].AqComp)
                    #print(self.grid[k, lastStep].IceComp)
                    for key in self.grid[k,lastStep].RockComp:
                        self.grid[k,i].RockComp[key]=self.grid[k,lastStep].RockComp[key]/sumC*self.grid[k,i].RockMass
        

    def updateBulk(self,time_step):
        '''
           Update the bulk properties of the model
        '''
        bot=0
        for i in range(0,self.nr):
            bot=self.grid[time_step,i].calc_thickness(bot)  # update the radii of the cells
            self.grid[time_step,i].calc_vol()               # update the volume of each cell
            self.grid[time_step,i].calc_Mass()              # update the mass of each cell
        self.calc_press(time_step)                          # update the pressure in each cell
    



    ##############################
    # Functions to be used by user
    ##############################
   
 
    def runModel(self):
        '''
           Run timestep function until end time or max time step is reached
           TODO: Write better doc string
        '''
        i=0
        while self.times[-1] < self.endTime and i<self.nt-1: #i in range(0,nt-1):
            self.timeStep(i)
            i=i+1
    

    def initialize_comp(self,name,initIceComp,initRockComp, initRIComp, initTemp, initRho, initK, initCp):
        '''
           Initializes the composition of the model
           TODO: Write better doc string
        '''
        self.name=name
        self.grid[0,0].IceComp=initIceComp
        self.grid[0,0].RockComp=initRockComp
        #self.grid[0,0].OrgComp=initOrgComp
        #self.grid[0,0].OrigOrgMass=initOrgMass
        self.grid[0,0].RIComp=initRIComp
        self.grid[0,0].Temp=initTemp
        self.grid[0,0].Dens=initRho
        self.grid[0,0].TCond=initK
        self.grid[0,0].Cp=initCp
        self.grid[0,0].TempStep=self.TempStep
        self.boundTemp=initTemp
        for i in range(1,self.nr):
            self.grid[0,i].copy_cell(self.grid[0,0],self.times[0])
            self.grid[0,i].ind=i
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


    ##### Plotting functions (ADD DOC STRINGS!)

    def plotAttribute(self,att,tit):
        matplotlib.rcParams.update({'font.size': 22})
        data = np.zeros((self.nt,self.nr))
        for i in range(self.nr):
            for j in range(self.nt):
                data[j][i]=getattr(self.grid[j][i], att)
        f=plt.figure(figsize=(10,6))

        xi=self.times/Decay.YR/10**6
        yi=self.radii[0:-1]/1000
        zi=np.transpose(data)

        CS1 = plt.contour( xi,yi,zi, 15, linewidths=2.0, colors='black', linestyles='dashed' )

        # labels on a subset of the major contour lines
        #labeled_levels = [ x for x in major_levels ]
        #clabels = plt.clabel( CS1, labeled_levels, fmt='%.0f', fontsize=12 )
        #for label in clabels:
        #    label.set_rotation(-90)

        CS3 = plt.contourf(xi,yi,zi, 500) #,vmin=vmin,vmax=vmax,cmap=JCmap) # Fix to use correct radii
        plt.colorbar()
        plt.xlabel('Time (Myr)')
        plt.ylabel('Radius (km)')
        plt.title(tit)
        return f
        

    def plotTemp(self,att,tit):
        matplotlib.rcParams.update({'font.size': 22})
        data = np.zeros((self.nt,self.nr))
        for i in range(self.nr):
            for j in range(self.nt):
                data[j][i]=getattr(self.grid[j][i], att)
        f=plt.figure(figsize=(10,6))
        #print(np.shape(self.radii[0:-1]))
        #print(np.shape(self.times))
        #print(np.shape(self.grid))
        vmin, vmax, JCmap = gmt_colormap()

        major_levels = [*range(300, 2300, 100)]
        #major_levels = [100,200,273] + major_levels
        major_levels = [100,200] + major_levels
        minor_levels = [ x for x in range(100, 1800, 25) if x not in major_levels ]
        xi=self.times/Decay.YR/10**6
        yi=self.radii[0:-1]/1000
        zi=np.transpose(data)

        CS1 = plt.contour( xi,yi,zi, levels=major_levels, linewidths=2.0, colors='black', linestyles='dashed' )

        # labels on a subset of the major contour lines
        labeled_levels = [ x for x in major_levels ]
        clabels = plt.clabel( CS1, labeled_levels, fmt='%.0f', fontsize=12 )
        for label in clabels:
            label.set_rotation(-90)

        CS2 = plt.contour( xi,yi,zi, levels=minor_levels, linewidths=1.0, colors='white', linestyles='dashed' )
        CS3 = plt.contourf(xi,yi,zi, 500,vmin=vmin,vmax=vmax,cmap=JCmap) # Fix to use correct radii
        plt.colorbar()
        #plt.contour(self.times/Decay.YR/10**6,self.radii[0:-1]/1000,np.transpose(data),15, colors='k')
        plt.xlabel('Time (Myr)')
        plt.ylabel('Radius (km)')
        plt.title(tit)
        return f
        
    def plotAttributeLine(self,att,tit,radius,unit=1):
        matplotlib.rcParams.update({'font.size': 22})
        data = np.zeros(self.nt)
        for j in range(self.nt):
            data[j]=getattr(self.grid[j][radius], att)
        f=plt.figure(figsize=(10,6))
        plt.plot(self.times/Decay.YR/10**6,data*unit,100) # Fix to use correct radii
        plt.xlabel('Time (Myr)')
        plt.ylabel(tit)
        #print(data)
        return f
    
    def plotDictAttribute(self,att,key,tit,unit=1):
        matplotlib.rcParams.update({'font.size': 22})
        data = np.zeros((self.nt,self.nr))
        for i in range(self.nr):
            for j in range(self.nt):
                tempItem=getattr(self.grid[j][i], att)
                if key in tempItem:
                    data[j][i]=tempItem[key]
        f=plt.figure(figsize=(10,6))
        plt.contourf(self.times/Decay.YR/10**6,self.radii[0:-1]/1000,np.transpose(data)*unit,100) # Fix to use correct radii
        plt.colorbar()
        plt.contour(self.times/Decay.YR/10**6,self.radii[0:-1]/1000,np.transpose(data)*unit,15, colors='k')
        plt.xlabel('Time (Myr)')
        plt.ylabel('Radius (km)')
        plt.title(tit)
        return f
        
    def plotDictAttributeMassScaled(self,att,key,tit,vmin,vmax,unit=1):
        matplotlib.rcParams.update({'font.size': 22})
        data = np.zeros((self.nt,self.nr))
        for i in range(self.nr):
            for j in range(self.nt):
                tempItem=getattr(self.grid[j][i], att)
                if key in tempItem:
                    data[j][i]=tempItem[key]/self.grid[j][i].Mass
        f=plt.figure(figsize=(10,6))
        plt.contourf(self.times/Decay.YR/10**6,self.radii[0:-1]/1000,np.transpose(data)*unit,100, vmin=vmin,vmax=vmax) # Fix to use correct radii
        plt.colorbar()
        plt.contour(self.times/Decay.YR/10**6,self.radii[0:-1]/1000,np.transpose(data)*unit,15, colors='k')
        plt.xlabel('Time (Myr)')
        plt.ylabel('Radius (km)')
        plt.title(tit)
        return f
    
    def plotDictAttributeLine(self,att,key,tit,radius,unit=1):
        matplotlib.rcParams.update({'font.size': 22})
        data = np.zeros(self.nt)
        for j in range(self.nt):
            tempItem=getattr(self.grid[j][radius], att)
            data[j]=0
            if key in tempItem:
                data[j]=tempItem[key]
        f=plt.figure(figsize=(10,6))
        plt.plot(self.times/Decay.YR/10**6,data*unit,100) # Fix to use correct radii
        plt.xlabel('Time (Myr)')
        plt.ylabel(tit)
        #print(data)
        return f,data

    def multiplotDictAttributeLine(self,att,keys,tit,radius,unit=1):
        matplotlib.rcParams.update({'font.size': 22})
        f = plt.figure(figsize=(10, 6))
        for key in keys:
            data = np.zeros(self.nt)
            for j in range(self.nt):
                tempItem=getattr(self.grid[j][radius], att)
                data[j]=0
                if key in tempItem:
                    data[j]=tempItem[key]
            plt.plot(self.times/Decay.YR/10**6,data*unit,100,label=key) # Fix to use correct radii
        plt.xlabel('Time (Myr)')
        plt.ylabel(tit)
        plt.legend()
        #print(data)
        return f


    def plotSimplePhases(self,radius,minorPhases="test"):
        pP = dict()

        if minorPhases=="test":
            minorPhases= ["Bugoff"]

        lab_dict = {"IOM": "IOM",
                    "IOMp": "IOM pyrolysates",
                    "other": "minor phases",
                    "gph": "Graphite",
                    "q": "Quartz",
                    "Opx": "Anhydrous silicates",
                    "Cpx": "Anhydrous silicates",
                    "acm": "Anhydrous silicates",
                    "pyr": "Sulfides",
                    "Po": "Sulfides",
                    "tro": "Sulfides",
                    "lot": "Sulfides",
                    "trov": "Sulfides",
                    "Amph": "Hydrous silicates",
                    "Chl": "Hydrous silicates",
                    "Tlc": "Hydrous silicates",
                    "Pu": "Hydrous silicates",
                    "Do": "Carbonates",
                    "Cc": "Carbonates",
                    "cc": "Carbonates",
                    "dol": "Carbonates",
                    "arag": "Carbonates",
                    "Pl": "Anhydrous silicates",
                    "ab": "Anhydrous silicates",
                    "hem": "Hematite",
                    "gth": "Goethite",
                    "Ol": "Anhydrous silicates",
                    "Sp": "Magnetite",
                    "Atg": "Hydrous silicates",
                    "glt": "Hydrous silicates",
                    "cen": "Anhydrous silicates",
                    "ne": "Anhydrous silicates",
                    "pren": "Anhydrous silicates",
                    "en": "Anhydrous silicates",
                    "liz": "Hydrous silicates",
                    "any": "Sulfates",
                    "naph": "Hydrous silicates",
                    "cor": "Al Oxides",
                    "Mag": "Carbonates",
                    "Gt": "Al Oxides",
                    "iron": "Iron",
                    "mt": "Magnetite"}

        col_dict = {"IOM": "#964B00",
                    "IOM pyrolysates": "#5C4033",
                    "minor phases": 'skyblue',
                    "Graphite": "#000000",
                    "Quartz": "#FFFFFF",
                    "Carbonates": "#808080",
                    "Anhydrous silicates": "#9ab973",
                    "Sulfides": "#E1C16E",
                    "Hydrous silicates": "#033220",
                    "Hematite": "#b7410e",
                    "Goethite": "#3b3b3b",
                    "Magnetite": "#5C4033",
                    "Al Oxides": "#9A2A2A",
                    "Sulfates": "#b0c4de",
                    "Iron": "#A9A9A9"}

        for j in range(self.nt):
            phases = self.grid[j][radius].RockPhases
            for i in range(len(phases)):
                if "_rs" in phases[i]:
                    phases[i] = phases[i].replace('_rs', '')
                if "_1" in phases[i]:
                    phases[i] = phases[i].replace('_1', '')
            for i in range(0, len(phases)):
                if not phases[i] == "Bulk":
                    labT = "other"
                    if phases[i] not in minorPhases and phases[i] not in lab_dict:
                        print(phases[i])
                    if phases[i] in lab_dict:
                        labT = lab_dict[phases[i]]
                    if labT in pP:
                        pP[labT][j] += self.grid[j][radius].RockPhaseDat[i]["wt%"]
                    else:
                        pP[labT] = np.zeros(self.nt)
                        pP[labT][j] = self.grid[j][radius].RockPhaseDat[i]["wt%"]
            if "IOM" in self.grid[j][radius].RockComp:
                if "IOM" in pP:
                    pP["IOM"][j] += self.grid[j][radius].RockComp["IOM"] / self.grid[j][radius].Mass * 100
                else:
                    pP["IOM"] = np.zeros(self.nt)
                    pP["IOM"][j] = self.grid[j][radius].RockComp["IOM"] / self.grid[j][radius].Mass * 100
        sumM = np.zeros(self.nt)
        for j in range(self.nt):
            for i in pP:
                if not i == "IOM":
                    sumM[j] += pP[i][j]
        for j in range(self.nt):
            if "IOM" in pP:
                IOMp = pP["IOM"][j]
            else:
                IOMp = 0
            for i in pP:
                if not i == "IOM":
                    pP[i][j] = pP[i][j] / sumM[j] * (100 - IOMp)

        # x_values = list(pP.keys())
        # y_values_list = list(pP.values())

        # print(pP.values())
        # print(pP.keys())

        color_map = []
        for i in pP:
            if i in col_dict:
                color_map.append(col_dict[i])
            else:
                color_map.append(col_dict['minor phases'])

        f = plt.figure(figsize=(12, 6))
        plt.stackplot(self.times / Decay.YR / 10 ** 6, pP.values(), labels=pP.keys(), colors=color_map)
        # plt.stackplot(x_values, y_values_list, labels=list(pP.keys()))
        plt.tight_layout(rect=[0, 0, 0.75, 1])
        # if LogScale==1:
        #    plt.xscale("log")
        plt.xlabel('Time (Myr)')
        plt.ylabel('Wt %')
        plt.title('Phase assemblage at %0.2f kms deep' % (self.radii[-1] - self.radii[radius]))
        plt.tick_params(labelright=True, right=True)
        plt.ylim([0, 100])
        plt.legend(bbox_to_anchor=(1.6, 1.0), loc='upper right',
                   ncol=1, fancybox=True, shadow=True)
        # plt.show()
        return f

   
    def plotPhaseAssemblage(self,radius, minorPhases="test"):
        pP=dict()
        
        if minorPhases=="test":
            minorPhases= ['F','Mica', 'Fsp', 'Sp', 'ky', 'Gt', 'trd', 'crst','Stlp','law','Pu','ank','Mica','Bio','Fsp','Sp']
        
        col_dict={"IOM":"#964B00",
                  "IOM pyrolysates":"#5C4033",
                  "minor phases":'skyblue',
                  "gph: graphite":"#000000",
                  "qtz: quartz":"#FFFFFF",
                  "cb: carbonates":"#808080",
                  "px: pyroxene":"#454B1B",
                  "py: pyrite":"#FFD700",
                  "po: pyrrhotite":"#E1C16E",
                  "amph: amphibole":"#033220",
                  "chl: chlorite":"#008080",
                  "tlc: talc":"#D3D3D3",
                  "pl: plagioclase":"#FF0000",
                  "hem: hematite":"#b7410e",
                  "gth: goethite":"#A9A9A9",
                  "mag: magnetite":"#5C4033",
                  "sp: spinel":"#E6E6FA",
                  "srp: serpentine":"#adb2a5",
                  "su: sulfates":"#b0c4de",
                  "phl: phlogopite":"#d53600",
                  "tro: troilite":"#C0C0C0",
                  "ol: olivine":"#9ab973",
                  "ens: enstatite":"#CC7722",
                  "cor: corundum":"#e0115f",
                  "Gt: garnet":"#9A2A2A",
                  "Iron":"#A9A9A9"}

        lab_dict={"IOM":"IOM",
                  "IOMp":"IOM pyrolysates",
                  "other":"minor phases",
                  "gph":"gph: graphite",
                  "q":"qtz: quartz",
                  "Opx":"px: pyroxene",
                  "Cpx":"px: pyroxene",
                  "acm":"px: pyroxene",
                  "pyr":"py: pyrite",
                  "Po":"po: pyrrhotite",
                  "tro":"tro: troilite",
                  "lot":"tro: troilite",
                  "trov":"po: pyrrhotite",
                  "Amph":"amph: amphibole",
                  "Chl":"chl: chlorite",
                  "Tlc":"tlc: talc",
                  "Do":"cb: carbonates",
                  "Cc":"cb: carbonates",
                  "cc": "cb: carbonates",
                  "dol":"cb: carbonates",
                  "arag":"cb: carbonates",
                  "Pl":"pl: plagioclase",
                  "ab":"pl: plagioclase",
                  "hem":"hem: hematite",
                  "gth": "gth: goethite",
                  "Ol":"ol: olivine",
                  "Sp":"sp: spinel",
                  "Atg":"srp: serpentine",
                  "glt":"srp: serpentine",
                  "cen":"ens: enstatite",
                  "pren": "ens: enstatite",
                  "en": "ens: enstatite",
                  "liz": "srp: serpentine",
                  "any":"su: sulfates",
                  "naph":"phl: phlogopite",
                  "cor":"cor: corundum",
                  "Mag":"cb: carbonates",#"mag: magnetite",
                  "Gt":"Gt: garnet",
                  "iron":"Iron"}


        for j in range(self.nt):
            phases=self.grid[j][radius].RockPhases
            for i in range(len(phases)):
                if "_rs" in phases[i]:
                    phases[i]=phases[i].replace('_rs','')
                if "_1" in phases[i]:
                    phases[i]=phases[i].replace('_1','')
            for i in range(0,len(phases)):
                if not phases[i]=="Bulk":
                    labT="other"
                    if phases[i] not in minorPhases and phases[i] not in lab_dict:
                        print(phases[i])
                    if phases[i] in lab_dict:
                        labT=lab_dict[phases[i]]
                    if labT in pP:
                        pP[labT][j]+=self.grid[j][radius].RockPhaseDat[i]["wt%"]
                    else:
                        pP[labT]=np.zeros(self.nt)
                        pP[labT][j]=self.grid[j][radius].RockPhaseDat[i]["wt%"]
            if "IOM" in self.grid[j][radius].RockComp:
                if "IOM" in pP: 
                    pP["IOM"][j]+=self.grid[j][radius].RockComp["IOM"]/self.grid[j][radius].Mass*100
                else:
                    pP["IOM"]=np.zeros(self.nt)
                    pP["IOM"][j]=self.grid[j][radius].RockComp["IOM"]/self.grid[j][radius].Mass*100
        sumM=np.zeros(self.nt)
        for j in range(self.nt):
            for i in pP:
                if not i=="IOM":
                    sumM[j]+=pP[i][j]
        for j in range(self.nt):
            if "IOM" in pP:
                IOMp=pP["IOM"][j]
            else:
                IOMp=0
            for i in pP:
                if not i=="IOM":
                    pP[i][j]=pP[i][j]/sumM[j]*(100-IOMp)

                
        #x_values = list(pP.keys())
        #y_values_list = list(pP.values())
   
        #print(pP.values())
        #print(pP.keys())



        color_map=[]
        for i in pP:
            if i in col_dict:
                color_map.append(col_dict[i])
            else:
                color_map.append(col_dict['minor phases'])
       

        f=plt.figure(figsize=(12,6))
        plt.stackplot(self.times/Decay.YR/10**6,pP.values(),labels=pP.keys(),colors = color_map) 
        #plt.stackplot(x_values, y_values_list, labels=list(pP.keys()))
        plt.tight_layout(rect=[0, 0, 0.75, 1])
        #if LogScale==1:
        #    plt.xscale("log")
        plt.xlabel('Time (Myr)')
        plt.ylabel('Wt %')
        plt.title('Phase assemblage at %0.2f kms deep'%(self.radii[-1]-self.radii[radius]))
        plt.tick_params(labelright=True, right=True)
        plt.ylim([0,100])
        plt.legend(bbox_to_anchor=(1.6, 1.0), loc='upper right',
                   ncol=1, fancybox=True, shadow=True)
        #plt.show()
        return f



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

def molMass2elMass(molDic):
    elDic={}
    for i in molDic:
        f = Formula(i)
        df = f.composition().dataframe()
        compDict=f.composition().dataframe().to_dict()["Fraction"]
        for j in compDict:
            if j in elDic:
                elDic[j]+=compDict[j]*molDic[i]
            else:
                elDic[j]=compDict[j]*molDic[i]
    return elDic


######### To be omitted

#def callAqFreeze(waterDict, iceDict, press, temp):
#    heat = 0
#    frozen=0
#    IceMeltHeat=334*1000 #J/kg
#    if temp<273.15:
#        if 'H2O' in waterDict:
#            iceDict['H2O']=waterDict['H2O']
#            heat=waterDict["H2O"]*IceMeltHeat
#            waterDict['H2O']=0
#            frozen=1
#    return heat,frozen
    
def callIceEquil(IceDict, waterDict, press, temp):
    heat = 0
    IceMeltHeat=334*1000 #J/kg
    if temp>273.15:
        if 'H2O' in IceDict:
            for i in IceDict:
                waterDict[i]=IceDict[i]
                IceDict[i]=0
            heat=waterDict["H2O"]*IceMeltHeat
    return heat
        














