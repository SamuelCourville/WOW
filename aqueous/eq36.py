import fileinput
import sys
import numpy as np
import os
import decimal
import matplotlib.pyplot as plt


eq36_dir="/Users/samuelcourville/Documents/JPL/combinedModel/aqueous/eq3_6/bin/"
eq3_ex=eq36_dir+"eq3nr"
eq6_ex=eq36_dir+"eq6"

eq36_db_dir="/Users/samuelcourville/Documents/JPL/combinedModel/aqueous/eq3_6/db/alphja39DBs/"

eqFileLoc="/Users/samuelcourville/Documents/JPL/combinedModel/aqueous/eq36_files/"
eq6file="/Users/samuelcourville/Documents/JPL/combinedModel/aqueous/eq36_files/ariel.6i"
eq3file="/Users/samuelcourville/Documents/JPL/combinedModel/aqueous/eq36_files/ariel.3i"
eq6ofile="/Users/samuelcourville/Documents/JPL/combinedModel/aqueous/eq36_files/ariel.6o"
eq3pfile="/Users/samuelcourville/Documents/JPL/combinedModel/aqueous/eq36_files/ariel.3p"


ctx = decimal.Context()
ctx.prec = 10

def updateTemp6i(temp):
    file6i = eq6file
    string = "|  [x] ( 0) Constant temperature:                                              |"
    scientific_notation = "{:.5e}".format(temp)
    stringReplace = "|             Value (C)         | "+scientific_notation+"| (tempcb)                        |\n"
    return replaceLine(file6i, string, stringReplace, 1)

def updateTemp3i(temp):
    file3i = eq3file
    string = "|Temperature (C)         |"
    scientific_notation = "{:.5e}".format(temp)
    stringReplace = "|Temperature (C)         | "+scientific_notation+"| (tempc)                                |\n"
    return replaceLine(file3i, string, stringReplace, 0)

def updatePress3i(P):
    file3i = eq3file
    string = "( 2) Value (bars)"
    scientific_notation = "{:.5e}".format(P)
    stringReplace = "|  [x] ( 2) Value (bars) | "+scientific_notation+"| (press)                                |\n"
    return replaceLine(file3i, string, stringReplace, 0)

def updatePress6i(P):
    file6i = eq6file
    string = "|  [x] ( 2) Constant pressure:                                                 |"
    scientific_notation = "{:.5e}".format(P)
    stringReplace = "|             Value (bars)      | "+scientific_notation+"| (pressb)                        |\n"
    return replaceLine(file6i, string, stringReplace, 1)


def updateReactant3i(reactant, molality):
    file3i = eq3file
    string = "|"+reactant
    scientific_notation = "{:.5e}".format(molality)
    stringReplace = "|"+reactant.ljust(48," ")+"| "+scientific_notation+"|Molality        |\n"
    return replaceLine(file3i, string, stringReplace, 0)

def updateReactant(reactant, moles):
    file6i = eq6file
    string = reactant
    scientific_notation = "{:.5e}".format(moles)
    stringReplace = "|->|Amount remaining (moles) | "+scientific_notation+"| (morr(n))                          |\n"
    stringReplace2 = "|--->|dXi(n)/dXi (mol/mol)      | "+scientific_notation+"| (rkb(1,1,n))                    |\n"
    replaceLine(file6i, string, stringReplace2, 22)
    return replaceLine(file6i, string, stringReplace, 6)

def replaceLine(file, lineStart,lineReplace,endCount):
    count = 0
    countStart=0
    out="not found"
    finput = fileinput.input(file, inplace=True)
    for line in finput:
        if count>0:
            if count == endCount:
                line = lineReplace
                sys.stdout.write(line)
                out="found"
            else:
                sys.stdout.write(line)
            count=count+1
        elif lineStart in line and "Couple (aux. sp.)" not in line:
            if endCount==0:
                line = lineReplace
                sys.stdout.write(line)
                out="found"
            else:
                count=1
                countStart=1
                sys.stdout.write(line)
        else:
            sys.stdout.write(line)
    fileinput.close()
    return(out)
           
    
def copyPickup(file6i,file3p):
    lineStart = "* Start of the bottom half of the input file                                   *"
    flag3p = 0
    finput = fileinput.input(file6i, inplace=True)
    for line in finput:
        if lineStart in line:
            flag3p = 1
        elif flag3p == 0:
            sys.stdout.write(line)
        else:
            pass
    fileinput.close()
    f1 = open(file6i, 'a+')
    f2 = open(file3p, 'r')
 
    # appending the contents of the second file to the first file
    f1.write(f2.read())
 
    # closing the files
    f1.close()
    f2.close()
    return

def extractVal(searchFor, retCurr):
    if os.path.exists("tab"):
        fp = open("tab")
    else:
        print("no tab!")
        return "not found"
    now=0
    previousline="nope"
    for line in fp:
        if now==1:
            return line
        if searchFor in line:
            if retCurr == 0:
                return previousline
            else:
                now=1
        previousline=line
    return("not found")

def numberOfSteps():
    count=0
    Xi = np.array([0])
    printActive=0
    try:
        with open(eq6ofile) as f:
            for line in f:
                if "Log Xi=" in line and printActive:
                    spl=line.split()
                    Xi = np.append(Xi,10**np.float64(spl[2]))
                    printActive=0
                if "--- Major Species by Contribution to Aqueous Mass Balances ---" in line:
                    count=count+1
                    printActive=1

                    
    except FileNotFoundError:
        print('sucks to suck')
    return count,Xi

def extractSpecie(spec,N):
    mSpec = np.zeros(N)
    count=-1
    try:
        with open(eq6ofile) as f:
            printActive=0
            for line in f:
                if printActive:
                    if spec in line:
                        splitLine = line.split()
                        mSpec[count]=np.float64(splitLine[1])
                        printActive=0
                    
                if "--- Major Species by Contribution to Aqueous Mass Balances ---" in line:
                    printActive=0
                if "--- Distribution of Aqueous Solute Species ---" in line:
                    printActive=1
                    count=count+1

    except FileNotFoundError:
        print('sucks to suck')
    return mSpec
        
def extractProduct(N):
    mSolids = dict()
    countdown=0
    count=-1
    try:
        with open(eq6ofile) as f:
            for line in f:
                if countdown>1:
                    countdown = countdown-1
                if countdown==1:
                    splitLine = line.split()
                    if np.size(splitLine)>0:
                        if splitLine[0] in mSolids:
                            mSolids[splitLine[0]][count]=np.float64(splitLine[1])
                        else:
                            mSolids[splitLine[0]]=np.zeros(N)*np.nan
                            mSolids[splitLine[0]][count]=np.float64(splitLine[1])
                    else:
                        countdown=0
                if "--- Grand Summary of Solid Phases (ES + PRS + Reactants) ---" in line:
                    countdown=5
                    count = count+1

    except FileNotFoundError:
        print('sucks to suck')
    return mSolids

def extractpH():
    mSpec = np.array([])
    off=1
    try:
        with open(eq6ofile) as f:
            for line in f:
                if "NBS pH scale" in line and off==1:
                    off=0
                elif "NBS pH scale" in line and off==0:
                    splitLine = line.split()
                    mSpec=np.append(mSpec,np.float64(splitLine[3]))
                    off=1
            
    except FileNotFoundError:
        print('sucks to suck')
    return mSpec
    
def moveFile(file,fileLoc,newLoc):
    status=os.system("mv "+fileLoc+file+" "+newLoc)
        
def runModel(tempK,Ppascal):
    P = Ppascal*10**-5 # convert to bars
    temp=tempK-273.15
    updateTemp3i(temp)
    updateTemp6i(temp)
    updatePress3i(P)
    updatePress6i(P)

    db500="500.d1"
    db1kb="1kb.d1"
    db2kb="2kb.d1"
    db5kb="5kb.d1"
    pwd="/Users/samuelcourville/Documents/JPL/combinedModel/"
    if P<500:
        status = os.system(eq3_ex +" "+eq36_db_dir+db500+" "+eq3file+" > /dev/null")
        moveFile("ariel.3p",pwd,eqFileLoc)
        copyPickup(eq6file,eq3pfile)
        status = os.system(eq6_ex +" "+eq36_db_dir+db500+" "+eq6file+" > /dev/null")
        moveFile("ariel.6o",pwd,eqFileLoc)
    elif P<1500:
        status = os.system(eq3_ex +" "+eq36_db_dir+db1kb+" "+eq3file+" > /dev/null")
        moveFile("ariel.3p",pwd,eqFileLoc)
        copyPickup(eq6file,eq3pfile)
        status = os.system(eq6_ex +" "+eq36_db_dir+db1kb+" "+eq6file+" > /dev/null")
        moveFile("ariel.6o",pwd,eqFileLoc)
    elif P<2500:
        status = os.system(eq3_ex +" "+eq36_db_dir+db2kb+" "+eq3file+" > /dev/null")
        moveFile("ariel.3p",pwd,eqFileLoc)
        copyPickup(eq6file,eq3pfile)
        status = os.system(eq6_ex +" "+eq36_db_dir+db2kb+" "+eq6file+" > /dev/null")
        moveFile("ariel.6o",pwd,eqFileLoc)
    else:
        status = os.system(eq3_ex +" "+eq36_db_dir+db5kb+" "+eq3file+" > /dev/null")
        moveFile("ariel.3p",pwd,eqFileLoc)
        copyPickup(eq6file,eq3pfile)
        status = os.system(eq6_ex +" "+eq36_db_dir+db5kb+" "+eq6file+" > /dev/null")
        moveFile("ariel.6o",pwd,eqFileLoc)
    
    #print("exit: " + str(status))
    return status

def updateAQ(H2,Cl):
    updateReactant3i("H2(aq)", H2)
    #updateReactant3i("HCO3-", HCO3)
    #updateReactant3i("NH3(aq)", NH3)
    updateReactant3i("Cl-", Cl)
    #updateReactant3i("Methanol(aq)", M)
    #updateReactant3i("Formaldehyde(aq)", F)

def updateRock(Q,S,C,K,N,L,B,F,P,Mn,CO2,CO,NH3,M,Fm,CH4):
    updateReactant("Quartz", Q)
    updateReactant("H2S(aq)", S)
    updateReactant("Corundum", C)
    updateReactant("K2O", K)
    updateReactant("Na2O", N)
    updateReactant("Lime", L)
    updateReactant("Bunsenite", B)
    updateReactant("FeO", F)
    updateReactant("Periclase", P)
    updateReactant("Manganosite", Mn)
    updateReactant("CO2(aq)", CO2)
    updateReactant("CO(aq)", CO)
    updateReactant("NH3(aq)", NH3)
    updateReactant("Methanol(aq)", M)
    updateReactant("Formaldehyde(aq)", Fm)
    updateReactant("Methane(aq)", CH4)

def extractAq(N):
    mAq = dict()
    countdown=0
    count=-1
    try:
        with open(eq6ofile) as f:
            for line in f:
                if countdown>1:
                    countdown = countdown-1
                if countdown==1:
                    splitLine = line.split()
                    if np.size(splitLine)>0:
                        if splitLine[0] in mAq:
                            mAq[splitLine[0]][count]=np.float64(splitLine[1])
                        else:
                            mAq[splitLine[0]]=np.zeros(N)*np.nan
                            mAq[splitLine[0]][count]=np.float64(splitLine[1])
                    else:
                        countdown=0
                if "--- Distribution of Aqueous Solute Species ---" in line:
                    countdown=5
                    count = count+1

    except FileNotFoundError:
        print('sucks to suck')
    return mAq

def extractFug(N):
    fugs = dict()
    countdown=0
    count=-1
    try:
        with open(eq6ofile) as f:
            for line in f:
                if countdown>1:
                    countdown = countdown-1
                if countdown==1:
                    splitLine = line.split()
                    if np.size(splitLine)>0:
                        if splitLine[0] in fugs:
                            try:
                                fugs[splitLine[0]][count]=np.float64(splitLine[2])
                            except ValueError:
                                fugs[splitLine[0]][count]=0
                        else:
                            fugs[splitLine[0]]=np.zeros(N)*np.nan
                            try:  
                                fugs[splitLine[0]][count]=np.float64(splitLine[2])
                            except ValueError: 
                                fugs[splitLine[0]][count]=0
                    else:
                        countdown=0
                if "--- Fugacities ---" in line:
                    countdown=5
                    count = count+1

    except FileNotFoundError:
        print('sucks to suck')
    return fugs

def extractEnd(diction):
    ends = dict()
    for i in diction:
        ends[i]=diction[i][-1]
    return ends

def float_to_str(f):
    """
    Convert the given float to a string,
    without resorting to scientific notation
    """
    if np.isnan(f):
        f=0.0
    d1 = ctx.create_decimal(repr(f))
    return format(d1, 'f')
