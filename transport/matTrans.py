# imports
import numpy as np
from aqStep import *

# Class with functions to handle moving fluids
class matTrans:

    def check_grid_mass_rock(grid,rockNew):
        '''
            Checks whether mass balance was conserved across rock cells after a mass transfer step.
                grid - array of grid cells
                rockNew - dictionary of rock compositions for each new cell
        '''
        mOld=0
        mNew=0
        for i in range(0,len(grid)):
            mOld+=grid[i].getRockMass()
            mNew+=matTrans.sumDictComps(rockNew[i])
        if abs(mOld-mNew)>1000000:
            print("rock imbalance")
            print(mOld)
            print(mNew)
        return False

    def check_grid_mass_water(grid,waterNew):
        '''
            Checks whether mass balance was conserved across fluid cells after a mass transfer step.
                grid - array of grid cells
                waterNew - dictionary of water compositions for each new cell
        '''
        mOld=0
        mNew=0
        for i in range(0,len(grid)):
            mOld+=grid[i].getAqMass()
            mNew+=matTrans.sumDictComps(waterNew[i])
        if abs(mOld-mNew)>1000000:
            print("water imbalance")
            print(mOld)
            print(mNew)
        return False

    def differentiate(grid,tC,bC):
        '''
            Moves fluid above rock. If Fluid is found in a cell beneath a rock cell, it moves that fluid upward.
            Unless that cell has porosity, in which case, this function only moves fluid if more than what would fit in
            pore space.
                grid - array of grid cells.
                tC - top most cell to consider for differentiation (do I still need this?)
                bC - bottom most cell to consider for differentiation (do I still need this?)
        '''
        nC=tC-bC

        # Remove all water from Ceres.
        aqAll = matTrans.grabWater(grid, nC)

        # Reconstruct Ceres with just rocks. Preserving order.
        rockCopy, phaseDatCopy, phasesCopy = matTrans.squishRocks(grid,nC)
        matTrans.check_grid_mass_rock(grid, rockCopy)

        # Distribute water in remaining cells.
        aqCopy= matTrans.distributeWater(grid,aqAll,rockCopy,nC)
        matTrans.check_grid_mass_water(grid, aqCopy)
        #orgCopy = matTrans.checkOrgs(grid,nC)  # Organics not worling right now
        return aqCopy,rockCopy,phaseDatCopy,phasesCopy

    def AddImpurities(aqcopy,impurities):
        '''
            Used to add elements to an aq dictionary
            When ice freezes in th ice shell, the remaining non-water elements are added to the fluid cell beneath.
                aqcopy - dictionary of aq elements to add impurities to
                impurities - Dictionary of elements remaining after ice shell freezing to add to fluid cell
        '''

        M=0
        for i in impurities:
            M+=impurities[i]


        for i in impurities:
            if i in aqcopy:
                aqcopy[i]+=impurities[i]
            else:
                aqcopy[i]=impurities[i]

        impurities={}
        return aqcopy,impurities



    def checkOrgs(grid,nC):
        '''
            Copies organic composition to new array
        '''
        orgCopy=[{}]*nC
        for i in range(0,nC):
             if not grid[i].getRockMass()==0:
                 orgCopy[i]=grid[i].OrgComp
        return orgCopy

    def grabWater(grid,nC):
        '''
            Loop through all grid cells and collect the total mass of fluid elements into one new dictionary
        '''
        aqAll=dict()
        for i in range(0,nC):
            for j in grid[i].AqComp:
                if j in aqAll:
                    aqAll[j]+=grid[i].AqComp[j]
                else:
                    aqAll[j]=grid[i].AqComp[j]
        return aqAll

    def distributeWater(grid,aqAll,rockCopy,nC):
        '''
            Look for cells that don't have enough mass, and fill the gap with water
                aqAll - all the fluid mass
                rockCopy - grid cells with just rocks
                nC - number of cells
        '''
        nR = len(grid)
        aqCopy=[{}] * nR
        for i in range(0,nC):
            aqCopy[i]=dict()
        weights = dict()
        sumAqs=0
        for i in aqAll:
            sumAqs+=aqAll[i]
        if sumAqs == 0:
            return aqCopy
        for i in aqAll:
            weights[i]=aqAll[i]/sumAqs

        for i in range(0,nC):
            rockM=0
            for j in rockCopy[i]:
                rockM+=rockCopy[i][j]
            cellM=grid[i].Mass-grid[i].getIceMass()
            if (cellM-rockM)>10000: # Why do I need a buffer here? Why is numerical precision a problem?
                addM=cellM-rockM
                for k in aqAll:
                    aqCopy[i][k]=weights[k]*addM
                    aqAll[k]=aqAll[k]-weights[k]*addM
        return aqCopy

    # Function to compress rocks after removing fluid. This is a mess
    def squishRocks(grid,nC):
        '''
            After extracting all the fluid, compress all rocks.
                TODO - Think about how porosity is handled
        '''
        nR = len(grid)          # number of cells
        rockCopy=[{}] * nR      # array of rock elemental abundances
        rockWeights=[{}] * nR   # dictionary of element mass fractions at each cell
        phaseDatCopy=[{}] * nR  # rock phase data
        phasesCopy=[[]] * nR    # rock phase names
        rockM=[0] * nR          # rock masses, total mass of rock in each cell
        rockM3=[0] * nR         # rock masses

        # loop to collect rock mass data
        for i in range(nC):     # loop over cells
            rockCopy[i]=dict()
            rockWeights[i]=dict()
            phaseDatCopy[i]=dict()
            phasesCopy[i]=[]
            for j in grid[i].RockComp:
                # make sure cell isn't undifferentiated.
                if not grid[i].Celltype==0: # UNDIFF
                    rockM[i]+=grid[i].RockComp[j]      # add cell rock element mass to total rock mass for cell
                    rockCopy[i][j]=grid[i].RockComp[j] # add element mass to rock mass dictionary copy.
            for j in grid[i].RockComp:
                if rockM[i]>0:
                    rockWeights[i][j]=grid[i].RockComp[j]/rockM[i] # get element mass fraction
                else:
                    rockWeights[i][j]=0 # zero if no rock
            rockM3[i]=rockM[i]          # copy rock mass list


        # Loop to redistribute rocks
        for i in range(nC):
            pf=1.0
            if grid[i].Celltype==3: # rock
                if grid[i].getRockMass()/grid[i].Mass>(1-grid[i].Porosity):
                    pf = grid[i].getRockMass()/grid[i].Mass
                else:
                    pf = (1-grid[i].Porosity)
            cellM=grid[i].Mass*pf #-grid[i].getIceMass() # cell mass that needs filled
            missM=cellM-rockM[i]  # missing mass
            if missM < -1e6:
                print(missM)
                print(error)
            if missM < 0:
                rockM[i]+=missM
                missM=0

            for j in range(i+1,nC):
                if rockM[j]>=missM:
                    for k in rockCopy[j]:
                        rockCopy[i][k]+=missM*rockWeights[j][k]
                        rockCopy[j][k]-=missM*rockWeights[j][k]
                        rockM[j]-=missM*rockWeights[j][k]

                    break
                if rockM[j]<missM:
                    for k in rockCopy[j]:
                        rockCopy[i][k]+=rockCopy[j][k]
                        rockCopy[j][k]=0
                    missM=missM-rockM[j]
                    rockM[j]=0
        rockM2=[0] * nR
        tempRI={}
        for i in range(nC):
            if grid[i].Celltype == 0:  # UNDIFF
                for j in grid[i].RIComp:
                    tempRI[j]=grid[i].RIComp[j]
            for j in rockCopy[i]:
                rockM2[i]+=rockCopy[i][j]
            if rockM2[i]>0:
                phaseDatCopy[i]=grid[i].RockPhaseDat.copy()
                phasesCopy[i]=grid[i].RockPhases.copy()

            if rockM3[i]>0:
                RIfact=rockM2[i]/rockM3[i]
                for j in grid[i].RIComp:
                    grid[i].RIComp[j]*=RIfact
            elif rockM3[i]==0 and rockM2[i]>0:
                for j in grid[i].RIComp:
                    grid[i].RIComp[j]=grid[i-1].RIComp[j]
            else:
                for j in grid[i].RIComp:
                    grid[i].RIComp[j]=0

            if grid[i].Celltype==0: # UNDIFF
                for k in grid[i].RockComp:
                    rockCopy[i][k]=grid[i].RockComp[k]
                for j in tempRI:
                    grid[i].RIComp[j]=tempRI[j]

            grid[i].reclassify()

        return rockCopy,phaseDatCopy,phasesCopy



    def differentiate_OLD(grid,topCell,bottom):
        '''
            obsolete
        '''
        mr = np.zeros(topCell-bottom)
        mw = np.zeros(topCell-bottom)
        mC = np.zeros(topCell-bottom)
        w = [{} for sub in range(bottom,topCell)]
        r = [{} for sub in range(bottom,topCell)]
        for i in range(bottom,topCell):
            mr[i-bottom] = grid[i].getRockMass()
            mw[i-bottom] = grid[i].getWaterMass()
            mC[i-bottom] = grid[i].Mass
            w[i-bottom] = grid[i].removeWaterMass()
            r[i-bottom] = grid[i].removeRockMass()
            wi=0
            ri=0
        print('start')
        print(bottom)
        print(topCell)
        for i in range(bottom,topCell):
            wcomp, rcomp, wi, ri=matTrans.grabMass(grid[i].Mass,mr,mw,w,r,wi,ri)
            print(i)
            print(wcomp)
            print(rcomp)
            grid[i].replaceMass(wcomp,rcomp)
            print(sum(grid[i].AqComp.values())+sum(grid[i].RockComp.values()))
            print(grid[i].Mass)
            grid[i].reclassify()
        print(bob)

    def extractPoreFluid(gridC):
        '''
            Looks through rock phases predicted by Perple_X and finds fluid phases
            Removes fluid phases from rock dictionary
        '''
        ex={}#{"H":0,"C":0,"Mg":0,"Al":0,"Si":0,"S":0,"Ca":0,"Fe":0,"O":0,"Na":0,"N":0}
        bv=0

        for fluid in ["F_rs","F_1_rs","F_2_rs","F_3_rs"]:
            if fluid in gridC.RockPhases:
                indP = np.where(gridC.RockPhases==fluid)
                ind = indP[0][0] # What is this garbage?
                for i in ["H","C","Mg","Al","Si","S","Ca","Fe","O","Na","N"]:
                    if i not in gridC.AqComp:
                        gridC.AqComp[i]=0

                    aToAdd=gridC.RockPhaseDat[ind][i]/100*gridC.RockPhaseDat[ind]['wt%']/100*gridC.Mass
                    if aToAdd>gridC.RockComp[i]:
                        aToAdd=gridC.RockComp[i]

                    gridC.AqComp[i]+=aToAdd
                    gridC.RockComp[i]-=aToAdd

                    if i in ex:
                        ex[i]+=aToAdd
                    else:
                        ex[i]=aToAdd
                    if gridC.RockComp[i]<0:
                        print(i)
                        print(gridC.RockComp[i])
                        print("negative rock mass in extract pore fluid")
                gridC.RockPhases=np.delete(gridC.RockPhases,ind)
                del gridC.RockPhaseDat[ind]

        return bv, ex

    def removeRockEl(gcell,ind):
        '''
            Remove element from rock dictionary. Obsolete?
        '''
        Els = ["H","C","Mg","Al","Si","S","Ca","Fe","O","Na","N"]
        for i in Els:
            removeI = gcell.RockPhaseDat[ind][i]/100*gcell.RockPhaseDat[ind]['wt%']/100*gcell.getRockMass()
            gcell.RockComp[i]=gcell.RockComp[i]-removeI

    def overturn():
        '''
            Potential future function to handle density stratification or convection in rock or ice shell
        '''
        return 1+1

    def grabMass(cellMass,Rmass,Wmass,water,rock,wi,ri):
        '''
            Obsolete?
        '''
        tempMass=0
        wToAdd={}
        rToAdd={}
        while tempMass<(cellMass-cellMass*0.9999) and ri<len(Rmass):
            if Rmass[ri]<=(cellMass-tempMass):
                tempMass=tempMass+Rmass[ri]
                Rmass[ri]=0
                rToAdd=matTrans.addDicts(rToAdd,rock[ri])
                rock[ri]={}
                ri=ri+1
            elif Rmass[ri]>(cellMass-tempMass):
                removeM=Rmass[ri]
                Rmass[ri]=Rmass[ri]-(cellMass-tempMass)
                removeM=removeM-Rmass[ri]
                tempMass=tempMass+removeM
                Rmasstemp,rock[ri]=matTrans.extractMass(rock[ri],removeM)
                rToAdd=matTrans.addDicts(rToAdd,Rmasstemp)
        while tempMass<(cellMass-cellMass*0.9999) and wi<len(Wmass):
            if Wmass[wi]<=(cellMass-tempMass):
                tempMass=tempMass+Wmass[wi]
                Wmass[wi]=0
                wToAdd=matTrans.addDicts(wToAdd,water[wi])
                water[wi]={}
                wi=wi+1
            elif Wmass[wi]>(cellMass-tempMass):
                removeW=Wmass[wi]
                Wmass[wi]=Wmass[wi]-(cellMass-tempMass)
                removeW=removeW-Wmass[wi]
                tempMass=tempMass+removeW
                Wmasstemp,water[wi]=matTrans.extractMass(water[wi],removeW)
                wToAdd=matTrans.addDicts(wToAdd,Wmasstemp)
                if cellMass-tempMass < 10000:
                    tempMass=cellMass        

        return wToAdd, rToAdd, wi, ri

    def extractMass(dict,MtoRemove):
        '''
            Remove mass from a mass dictionary
        '''
        sumMass = 0
        for i in dict:
            sumMass=sumMass+dict[i]
        ratio = MtoRemove/sumMass
        dict2=dict.copy()
        for i in dict2:
            dict2[i]=dict2[i]*ratio
            dict[i]=dict[i]-dict2[i]
        return dict2, dict

    def sumDictComps(dict):
        '''
            sum the elements of a dictionary together
        '''
        sumMs=0
        for i in dict:
            sumMs += dict[i]
        return sumMs
        
    def addDicts(x,y):
        '''
            Add elements of two dictionaries together
        '''
        return {k: x.get(k, 0) + y.get(k, 0) for k in set(x) | set(y)}

