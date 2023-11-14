import numpy as np

class matTrans:

    def differentiate(grid,tC,bC):
        nC=tC-bC
        
        rockCopy, phaseDatCopy, phasesCopy = matTrans.squishRocks(grid,nC)
        aqAll = matTrans.grabWater(grid,nC)
        aqCopy= matTrans.distributeWater(grid,aqAll,rockCopy,nC)
        return aqCopy,rockCopy,phaseDatCopy,phasesCopy

    def grabWater(grid,nC):
        aqAll=dict()
        for i in range(0,nC):
            for j in grid[i].AqComp:
                if j in aqAll:
                    aqAll[j]+=grid[i].AqComp[j]
                else:
                    aqAll[j]=grid[i].AqComp[j]
        return aqAll

    def distributeWater(grid,aqAll,rockCopy,nC):
        nR = len(grid)
        aqCopy=[{}] * nR
        for i in range(0,nC):
            aqCopy[i]=dict()
        weights = dict()
        sumAqs=0
        for i in aqAll:
            sumAqs+=aqAll[i]
        for i in aqAll:
            weights[i]=aqAll[i]/sumAqs

        for i in range(0,nC):
            rockM=0
            for j in rockCopy[i]:
                rockM+=rockCopy[i][j]
            cellM=grid[i].Mass
            if rockM<cellM:
                addM=cellM-rockM
                for k in aqAll:
                    aqCopy[i][k]=weights[k]*addM
                    aqAll[k]=aqAll[k]-weights[k]*addM
            #print(rockM)
            #print(cellM)
            #print(" ")
        #print(bob)
        return aqCopy

    def squishRocks(grid,nC):
        nR = len(grid)
        rockCopy=[{}] * nR
        rockWeights=[{}] * nR
        phaseDatCopy=[{}] * nR
        phasesCopy=[[]] * nR
        rockM=[0] * nR
        rockM3=[0] * nR
        for i in range(nC):
            rockCopy[i]=dict()
            rockWeights[i]=dict()
            phaseDatCopy[i]=dict()
            phasesCopy[i]=[]
            for j in grid[i].RockComp:
                rockM[i]+=grid[i].RockComp[j]
                rockCopy[i][j]=grid[i].RockComp[j]
            for j in grid[i].RockComp:
                if rockM[i]>0:
                    rockWeights[i][j]=grid[i].RockComp[j]/rockM[i]
                else:
                    rockWeights[i][j]=0
            rockM3[i]=rockM[i]
        for i in range(nC):
            cellM=grid[i].Mass
            missM=cellM-rockM[i]
            addRock=dict()
            for j in range(i+1,nC):
                if rockM[j]>=missM:
                    for k in rockCopy[j]:
                        rockCopy[i][k]+=missM*rockWeights[j][k]
                        rockCopy[j][k]-=missM*rockWeights[j][k]
                        rockM[j]-=missM*rockWeights[j][k]
                        #print(rockCopy[j][k])
                    #print(rockCopy[i]["H"])
                    #print("1")
                    #print("")
                    missM=0
                    break
                if rockM[j]<missM:
                    for k in rockCopy[j]:
                        rockCopy[i][k]+=rockCopy[j][k]
                        rockCopy[j][k]=0
                    missM=missM-rockM[j]
                    rockM[j]=0
                    #print(rockCopy[i]["H"])
                    #print("2")
            #print(rockCopy[i]["H"])
        #print(bob)
        rockM2=[0] * nR
        for i in range(nC):
            #print(grid[i].Mass)
            for j in rockCopy[i]:
                rockM2[i]+=rockCopy[i][j]
            #print(rockM2[i])
            #print(" ")
            if rockM2[i]>0:
                phaseDatCopy[i]=grid[i].RockPhaseDat.copy()
                phasesCopy[i]=grid[i].RockPhases.copy()
            #else:
            #    phasesCopy[i]=grid[i].RockPhases.copy()
            #    phaseDatCopy[i]=grid[i].RockPhaseDat.copy()

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
            
            grid[i].reclassify()
                
        return rockCopy,phaseDatCopy,phasesCopy



    def differentiate_OLD(grid,topCell,bottom):
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
        if "F_rs" in gridC.RockPhases:
            indP = np.where(gridC.RockPhases=='F_rs')
            ind = indP[0][0] # What is this garbage?
            #print(ind)
            #print(gridC.RockPhases)
            #print(len(gridC.RockPhaseDat))
            #print(gridC.RockPhaseDat)
            if 'H2O' in gridC.AqComp:
                gridC.AqComp['H2O']+=gridC.RockPhaseDat[ind]['wt%']/100*gridC.getRockMass()
            else:                
                gridC.AqComp['H2O']=gridC.RockPhaseDat[ind]['wt%']/100*gridC.getRockMass()
            matTrans.removeRockEl(gridC,ind)
            gridC.RockPhases=np.delete(gridC.RockPhases,ind)
            del gridC.RockPhaseDat[ind]
            #print(bob)
        if "F_1_rs" in gridC.RockPhases:
            indP = np.where(gridC.RockPhases=='F_1_rs')
            ind = indP[0][0] # What is this garbage?
            if 'H2O' in gridC.AqComp:
                gridC.AqComp['H2O']+=gridC.RockPhaseDat[ind]['wt%']/100*gridC.getRockMass()
            else:                
                gridC.AqComp['H2O']=gridC.RockPhaseDat[ind]['wt%']/100*gridC.getRockMass()
            matTrans.removeRockEl(gridC,ind)
            gridC.RockPhases=np.delete(gridC.RockPhases,ind)
            del gridC.RockPhaseDat[ind]
        if "F_2_rs" in gridC.RockPhases:
            print("found F_2, but no code to deal with it yet.")

    def removeRockEl(gcell,ind):
        Els = ["H","C","Mg","Al","Si","S","Ca","Fe","O","Na","N"]
        for i in Els:
            #print("")
            #print(i)
            #print(gcell.RockComp[i])
            removeI = gcell.RockPhaseDat[ind][i]/100*gcell.RockPhaseDat[ind]['wt%']/100*gcell.getRockMass()
            gcell.RockComp[i]=gcell.RockComp[i]-removeI
            #print(removeI)
            #print(gcell.RockComp[i])

    def overturn():
        return 1+1

    def grabMass(cellMass,Rmass,Wmass,water,rock,wi,ri):
        tempMass=0
        wToAdd={}
        rToAdd={}
        while tempMass<(cellMass-cellMass*0.9999) and ri<len(Rmass):
            if Rmass[ri]<=(cellMass-tempMass):
                #print('1111')
                #print(ri)
                #print(cellMass)
                #print(Rmass[ri])
                #print(rock[ri])
                tempMass=tempMass+Rmass[ri]
                Rmass[ri]=0
                rToAdd=matTrans.addDicts(rToAdd,rock[ri])
                rock[ri]={}
                ri=ri+1
                #print(ri)
                #print(Rmass[ri-1])
                #print(rock[ri-1])
            elif Rmass[ri]>(cellMass-tempMass):
                #print('2222')
                #print(ri)
                #print(cellMass)
                #print(Rmass[ri])
                #print(rock[ri])
                removeM=Rmass[ri]
                Rmass[ri]=Rmass[ri]-(cellMass-tempMass)
                removeM=removeM-Rmass[ri]
                tempMass=tempMass+removeM
                Rmasstemp,rock[ri]=matTrans.extractMass(rock[ri],removeM)
                rToAdd=matTrans.addDicts(rToAdd,Rmasstemp)
                #print(ri)
                #print(Rmass[ri])
                #print(rock[ri])
        while tempMass<(cellMass-cellMass*0.9999) and wi<len(Wmass):
            if Wmass[wi]<=(cellMass-tempMass):
                #print('3')
                tempMass=tempMass+Wmass[wi]
                Wmass[wi]=0
                wToAdd=matTrans.addDicts(wToAdd,water[wi])
                water[wi]={}
                wi=wi+1
            elif Wmass[wi]>(cellMass-tempMass):
                #print('4')
                #print(Wmass[wi])
                #print(cellMass-tempMass)
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
        sumMass = 0
        for i in dict:
            sumMass=sumMass+dict[i]
        ratio = MtoRemove/sumMass
        dict2=dict.copy()
        for i in dict2:
            dict2[i]=dict2[i]*ratio
            dict[i]=dict[i]-dict2[i]
        return dict2, dict
        
    def addDicts(x,y):
        return {k: x.get(k, 0) + y.get(k, 0) for k in set(x) | set(y)}

