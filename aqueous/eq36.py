import fileinput
import sys
import numpy as np
import os
import decimal
import platform
import matplotlib.pyplot as plt

import WOW

##### UPDATE THIS LINE
modDir=WOW.main_directory
eq36Dir=WOW.eq36_directory


# File locations
eq36_dir=eq36Dir+"bin/"
eq3_ex=eq36_dir+"eq3nr"
eq6_ex=eq36_dir+"eq6"
eq36_db_dir=eq36Dir+"db/alphja39DBs/"
eqFileLoc=modDir+"aqueous/eq36_files/"
eq6master=modDir+"aqueous/eq36_files/master.6i"
eq3master=modDir+"aqueous/eq36_files/master.3i"

ctx = decimal.Context()
ctx.prec = 10


def checkEQfileExist(name,ind):
    '''
        Make sure EQ6 input file exists before running
    '''
    return os.path.isfile(eqFileLoc+name+str(ind)+".6i")
    


def copy_text_file(original_file, copied_file):
    '''
        Copy a text file
    '''
    try:
        with open(original_file, 'r') as original:
            # Read the contents of the original file
            content = original.read()

            # Write the contents to the new file
            with open(copied_file, 'w') as copied:
                copied.write(content)

        #print(f"Contents of '{original_file}' successfully copied to '{copied_file}'.")
    except FileNotFoundError:
        print(f"Error: The file '{original_file}' does not exist.")
    except Exception as e:
        print(f"An error occurred: {e}")

def copyMasters(name,r):
    '''
        Copies the template EQ36 input files into new files specific for certain grid cells
    '''
    #print("t1")
    copy_text_file(eq6master, eqFileLoc+name+str(r)+'.6i')
    copy_text_file(eq3master, eqFileLoc+name+str(r)+'.3i')

def rename6p_to_6i(file6p,file6i):
    '''
        Change EQ6 pickup file to input file
    '''
    #print("t2")
    copy_text_file(file6p, file6i)


def updateTemp6i(temp,file6i):
    '''
        Update temperature in EQ6i file
    '''
    #file6i = eq6file
    string = "|  [x] ( 0) Constant temperature:                                              |"
    scientific_notation = "{:.5e}".format(temp)
    stringReplace = "|             Value (C)         | "+scientific_notation+"| (tempcb)                        |\n"
    return replaceLine(file6i, string, stringReplace, 1)

def updateTemp3i(temp,file3i):
    '''
        Update temperature in EQ3i file
    '''
    #file3i = eq3file
    string = "|Temperature (C)         |"
    scientific_notation = "{:.5e}".format(temp)
    stringReplace = "|Temperature (C)         | "+scientific_notation+"| (tempc)                                |\n"
    return replaceLine(file3i, string, stringReplace, 0)

def updatePress3i(P,file3i):
    '''
        Update pressure in EQ3i file
    '''
    #file3i = eq3file
    string = "( 2) Value (bars)"
    scientific_notation = "{:.5e}".format(P)
    stringReplace = "|  [x] ( 2) Value (bars) | "+scientific_notation+"| (press)                                |\n"
    return replaceLine(file3i, string, stringReplace, 0)

def updatePress6i(P,file6i):
    '''
        Update pressure in EQ6i file
    '''
    #file6i = eq6file
    string = "|  [x] ( 2) Constant pressure:                                                 |"
    scientific_notation = "{:.5e}".format(P)
    stringReplace = "|             Value (bars)      | "+scientific_notation+"| (pressb)                        |\n"
    return replaceLine(file6i, string, stringReplace, 1)


def updateReactant3i(reactant, molality, name):
    '''
        Update reactant in EQ3i file
    '''
    file3i = eqFileLoc+name+".3i"
    string = "|"+reactant
    scientific_notation = "{:.5e}".format(molality)
    stringReplace = "|"+reactant.ljust(48," ")+"| "+scientific_notation+"|Molality        |\n"
    return replaceLine(file3i, string, stringReplace, 0)

def updateReactant(reactant, moles,name):
    '''
        Update reactant in EQ6i file
    '''
    file6i = eqFileLoc+name+".6i"
    string = reactant
    scientific_notation = "{:.5e}".format(moles)
    stringReplace = "|->|Amount remaining (moles) | "+scientific_notation+"| (morr(n))                          |\n"
    stringReplace2 = "|--->|dXi(n)/dXi (mol/mol)      | "+scientific_notation+"| (rkb(1,1,n))                    |\n"
    replaceLine(file6i, string, stringReplace2, 22)
    return replaceLine(file6i, string, stringReplace, 6)

def replaceLine(file, lineStart,lineReplace,endCount):
    '''
        Replace line in a text file
    '''
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
    '''
        Copy pickup EQ3p file onto end of EQ6i file
    '''
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
    '''
        Obsolete
    '''
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

def numberOfSteps(name):
    '''
        Return number of steps taken in EQ6 calculation
    '''
    count=0
    Xi = np.array([0])
    printActive=0
    try:
        with open(eqFileLoc+name+'.6o') as f:
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

def extractSpecie(spec,N,name):
    '''
        Extract solute specie from EQ6 output
    '''
    mSpec = np.zeros(N)
    count=-1
    try:
        with open(eqFileLoc+name+'.6o') as f:
            printActive=0
            for line in f:
                if printActive:
                    if spec in line:
                        splitLine = line.split()
                        mSpec[count]=np.float64(splitLine[1])
                        #aSpec[count]=np.float64(splitLine[4])
                        printActive=0
                    
                if "--- Major Species by Contribution to Aqueous Mass Balances ---" in line:
                    printActive=0
                if "--- Distribution of Aqueous Solute Species ---" in line:
                    printActive=1
                    count=count+1

    except FileNotFoundError:
        print('sucks to suck')
    return mSpec
        
def extractProduct(N,name):
    '''
        Extract product mineral from EQ6 output
    '''
    mSolids = dict()
    countdown=0
    count=-1
    try:
        with open(eqFileLoc+name+'.6o') as f:
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

def extractpH(name):
    '''
        Extract pH from EQ6 output file
    '''
    mSpec = np.array([])
    off=1
    try:
        with open(eqFileLoc+name+'.6o') as f:
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
    '''
        Move file to new location
    '''
    if platform.system() == "Windows":
        status = os.system("move " + fileLoc + file + " " + newLoc)
    else:
        status=os.system("mv "+fileLoc+file+" "+newLoc)


def executeEQ6(tempK,Pa,name,ind):
    '''
        Execute EQ6 input file
    '''
    P = Pa*10**-5 # convert to bars
    temp=tempK-273.15
    
    eq6file=name+ind
    eq6ifile=eqFileLoc+name+ind+".6i"
    eq6pfile=eqFileLoc+name+ind+".6p"
    rename6p_to_6i(eq6pfile,eq6ifile)
    
    updateTemp6i(temp,eq6ifile)
    updatePress6i(P,eq6ifile)

    db500="500.d1"
    db1kb="1kb.d1"
    db2kb="2kb.d1"
    db5kb="5kb.d1"
    #pwd="/Users/samuelcourville/Documents/JPL/combinedModel/" # dangerous!
    pwd=modDir
    
    dbStr='5kb.d1'
    if P<500:
        dbStr='500.d1'
    elif P<1500:
        dbStr='1kb.d1'
    elif P<2500:
        dbStr='2kb.d1'

    if platform.system() == "Windows":
        status = os.system(eq6_ex +" "+eq36_db_dir+dbStr+" "+eq6ifile)
    else:
        status = os.system(eq6_ex + " " + eq36_db_dir + dbStr + " " + eq6ifile + " > /dev/null")
    moveFile(eq6file+".6o",pwd,eqFileLoc)
    moveFile(eq6file+".6tx",pwd,eqFileLoc)
    moveFile(eq6file+".6t",pwd,eqFileLoc)
    moveFile(eq6file+".6p",pwd,eqFileLoc)
    moveFile(eq6file+".6bb",pwd,eqFileLoc)
    moveFile(eq6file+".6ba",pwd,eqFileLoc)

    #print("exit: " + str(status))
    return status

 
def runModel(tempK,Ppascal,name,ind):
    '''
        Obsolete
    '''
    P = Ppascal*10**-5 # convert to bars
    temp=tempK-273.15
    
    eqfile=name+ind
    eq3ifile=eqFileLoc+name+ind+".3i"
    eq3pfile=eqFileLoc+name+ind+".3p"
    eq6ifile=eqFileLoc+name+ind+".6i"
    
    updateTemp3i(temp,eq3ifile)
    updateTemp6i(temp,eq6ifile)
    updatePress3i(P,eq3ifile)
    updatePress6i(P,eq6ifile)

    db500="500.d1"
    db1kb="1kb.d1"
    db2kb="2kb.d1"
    db5kb="5kb.d1"
    #pwd="/Users/samuelcourville/Documents/JPL/combinedModel/"
    #pwd="/Users/samuelcourville/Documents/JPL/combinedModel/aqueous/"
    pwd = modDir
    
    dbStr='5kb.d1'
    if P<500:
        dbStr='500.d1'
    elif P<1500:
        dbStr='1kb.d1'
    elif P<2500:
        dbStr='2kb.d1'

    if platform.system() == "Windows":
        status = os.system(eq3_ex +" "+eq36_db_dir+dbStr+" "+eq3ifile)
    else:
        status = os.system(eq3_ex + " " + eq36_db_dir + dbStr + " " + eq3ifile + " > /dev/null")
    moveFile(eqfile+".3p",pwd,eqFileLoc)
    moveFile(eqfile+".3o",pwd,eqFileLoc)
    copyPickup(eq6ifile,eq3pfile)
    if platform.system() == "Windows":
        status = os.system(eq6_ex +" "+eq36_db_dir+dbStr+" "+eq6ifile)
    else:
        status = os.system(eq6_ex + " " + eq36_db_dir + dbStr + " " + eq6ifile + " > /dev/null")
    moveFile(eqfile+".6o",pwd,eqFileLoc)
    moveFile(eqfile+".6tx",pwd,eqFileLoc)
    moveFile(eqfile+".6t",pwd,eqFileLoc)
    moveFile(eqfile+".6p",pwd,eqFileLoc)
    moveFile(eqfile+".6bb",pwd,eqFileLoc)
    moveFile(eqfile+".6ba",pwd,eqFileLoc)
    
    #print("exit: " + str(status))
    return status

def updateAQ(H2,Cl, dir):
    '''
        Update H2 and Cl in EQ3i file
    '''
    updateReactant3i("H2(aq)", H2, dir)
    updateReactant3i("Cl-", Cl, dir)

# MUST UPGRADE
def updateReacts(rs,name):
    '''
        Update reactants. Obsolete?
    '''
    for i in rs:
        if not i=="H2O":
            updateReactant(convert2EQ(i),rs[i],name)

# MUST UPGRADE!!!!!
def convert2EQ(name):
    '''
        Add (aq) to end of specie name
    '''
    return name+"(aq)"

def updateRock(Q,S,C,K,N,L,B,F,P,Mn,CO2,CO,NH3,M,Fm,CH4,name):
    '''
        Update specific reactants
        Obsolete
    '''
    updateReactant("Quartz", Q,name)
    updateReactant("H2S(aq)", S,name)
    updateReactant("Corundum", C,name)
    updateReactant("K2O", K,name)
    updateReactant("Na2O", N,name)
    updateReactant("Lime", L,name)
    updateReactant("Bunsenite", B,name)
    updateReactant("FeO", F,name)
    updateReactant("Periclase", P,name)
    updateReactant("Manganosite", Mn,name)
    updateReactant("CO2(aq)", CO2,name)
    updateReactant("CO(aq)", CO,name)
    updateReactant("NH3(aq)", NH3,name)
    updateReactant("Methanol(aq)", M,name)
    updateReactant("Formaldehyde(aq)", Fm,name)
    updateReactant("Methane(aq)", CH4,name)

def extractAq(N,name):
    '''
        Extract specie from EQ6 file
    '''
    mAq = dict()
    aAq = dict()
    countdown=0
    count=-1
    try:
        with open(eqFileLoc+name+'.6o') as f:
            for line in f:
                if countdown>1:
                    countdown = countdown-1
                if countdown==1:
                    splitLine = line.split()
                    if np.size(splitLine)>0:
                        if splitLine[0] in mAq:
                            mAq[splitLine[0]][count]=np.float64(splitLine[1])
                            aAq[splitLine[0]][count]=np.float64(splitLine[4])
                        else:
                            mAq[splitLine[0]]=np.zeros(N)*np.nan
                            mAq[splitLine[0]][count]=np.float64(splitLine[1])
                            aAq[splitLine[0]] = np.zeros(N) * np.nan
                            aAq[splitLine[0]][count] = np.float64(splitLine[4])
                    else:
                        countdown=0
                if "--- Distribution of Aqueous Solute Species ---" in line:
                    countdown=5
                    count = count+1

    except FileNotFoundError:
        print('sucks to suck')
    return mAq, aAq

def extractFug(N,name):
    '''
        Extract gas fugacity from EQ6 file
    '''
    fugs = dict()
    countdown=0
    count=-1
    try:
        with open(eqFileLoc+name+'.6o') as f:
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
    '''
        Extract last concentration value from an array of concentrations along EQ6 reaction output
    '''
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

def updateSpecialReactantElement(El,value,fname):
    '''
        Updates the elemental abundance of one element of an EQ6 special reactant
    '''
    filename=fname
    try:
        with open(filename, 'r') as file:
            lines = file.readlines()
        
        for i, line in enumerate(lines):
            if "|--->|"+El+" " in line and "(uesri(i,n), cesri(i,n))" in line:

                # Replace the numeric value with the new value
                if len(El)>1:
                    new_line = "|--->|"+El+"      |"+"{:0.16e}".format(value)+"| (uesri(i,n), cesri(i,n))                |\n"
                else:
                    new_line = "|--->|"+El+"       |"+"{:0.16e}".format(value)+"| (uesri(i,n), cesri(i,n))                |\n"
                # Update the line in the list of lines
                lines[i] = new_line

        with open(filename, 'w') as file:
            file.writelines(lines)
        #print(f"Successfully replaced the numeric value with {new_value}.")
    except Exception as e:
        print(f"Error: {e}")


# Obsoslete
def reweightSpecialReactant(values,massW):
    '''
        Obsolete
    '''
    newDict={}
    elD={"H":1,"C":12.01,"Mg":24.31,"Si":28.09,"Fe":55.85,"Ca":40.08,"Na":22.99,"Al":26.98,"K":39.10,"O":16,"S":32.07,"N":14.01}
    sumM=0
    sumKG=0
    for i in elD:
        if i in values:
            sumKG+=values[i]
            newDict[i]=values[i]/elD[i]
            sumM+=newDict[i]
    newM=0
    if sumM==0:
        return {},0.0
    for i in newDict:
        newDict[i]=newDict[i]/sumM
        newM+=newDict[i]*elD[i]
    newM=newM/1000
    moles=sumKG/massW/newM
    return newDict, moles
    

def updateSpecialReactant(values,ni):
    '''
        Updates the elemental abundances of an EQ6 special reactant
    '''
    name = eqFileLoc+ni+".6i"
    elementList=["H","C","S","N","Mg","Si","Fe","Ca","Na","Al","K","O"]
    #molVal, nMols=reweightSpecialReactant(values,massW) # DELETE
    for i in elementList:
        if i in values:
            updateSpecialReactantElement(i,values[i],name)
        else:
            updateSpecialReactantElement(i,0.0,name)
    updateSpecReactant("d0001",1.0,name)


def updateSpecReactant(reactant, moles,name):
    '''
        Appendage of above function
    '''
    file6i = name
    string = reactant
    scientific_notation = "{:.5e}".format(moles)
    stringReplace = "|->|Amount remaining (moles) | "+scientific_notation+"| (morr(n))                          |\n"
    stringReplace2 = "|--->|dXi(n)/dXi (mol/mol)      | "+scientific_notation+"| (rkb(1,1,n))                    |\n"
    replaceLine(file6i, string, stringReplace2, 46)
    return replaceLine(file6i, string, stringReplace, 6)

