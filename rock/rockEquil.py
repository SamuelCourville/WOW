import os
import rpy2.robjects as robjects
import numpy as np

#codeDir = "/Users/samuelcourville/Documents/JPL/combinedModel/rock/Rcrust/code/"
codeDir = "/Users/samuelcourville/Documents/JPL/Perplex/Rcrust/code/"
inputDir="/Users/samuelcourville/Documents/JPL/Perplex/Rcrust/Projects/WOW/Inputs/WOW.txt"
#inputDir="/Users/samuelcourville/Documents/JPL/combinedModel/rock/Rcrust/Projects/WOW/Inputs/WOW.txt"
mainScript="/Users/samuelcourville/Documents/JPL/Perplex/Rcrust/code/main.r"
#mainScript="/Users/samuelcourville/Documents/JPL/combinedModel/rock/Rcrust/code/main.r"
rockEquilDir="/Users/samuelcourville/Documents/JPL/combinedModel/rock"
fileRData = "/Users/samuelcourville/Documents/JPL/Perplex/Rcrust/Projects/WOW/WOW.RData"
#fileRData = "/Users/samuelcourville/Documents/JPL/combinedModel/rock/Rcrust/Projects/WOW/WOW.RData"

def rockEquil(Comp,P,T):
    fN=inputDir
    #fN="/Users/samuelcourville/Documents/JPL/combinedModel/rock/Rcrust/Projects/WOW/Inputs/WOW.txt"
    RCrustProjName="WOW"
    Pc = P/100000/1000 #convert to kbar
    sumComp = 0
    for key in Comp:    
        sumComp += Comp[key]
    for key in Comp:    
        Comp[key] *= 100/sumComp
    updateRCrustInput(fN,Pc,T,Comp)
    executeRCrust(RCrustProjName)
    ddicts, specs = extractRData()
    spec_list = list(specs)
    specAr = np.array(spec_list)
    success=0
    E=0
    if "Bulk_rs" in specs:
        success=1
        indx=specs.index("Bulk_rs")
        E = ddicts[indx]['Enthalpy (J/kg)']
        newCompDict=copy_dict_entries(ddicts[indx],Comp.keys())
        for key in newCompDict:    
            newCompDict[key] *= (1/100)
    else:
        newCompDict={}
        print("Perplex failed on input:")
        print(Comp)
    return newCompDict, ddicts, specAr, E, success

def copy_dict_entries(original_dict, keys_to_copy):
    new_dict = {}
    for key in keys_to_copy:
        if key in original_dict:
            new_dict[key] = original_dict[key]
    return new_dict



def executeRCrust(pN):
    #os.chdir("/Users/samuelcourville/Documents/JPL/combinedModel/rock/Rcrust/code/")
    os.chdir(codeDir)
    status = os.system("Rscript "+mainScript+" "+pN+" >/dev/null")
    os.chdir(rockEquilDir)

def updateRCrustInput(fN,P,T,Comp):
    #line 16 and 17
    line16="x_n<-1"
    line17="y_n<-1"
    replace_line_in_file(fN, 16, line16, fN)
    replace_line_in_file(fN, 17, line17, fN)

    #line 25
    line25="pt_definitions<-list(\"{1;1}_{1;1}\"=c(\""+str(P)+"+0*y_i\",\""+str(T)+"+0*x_i\"))"
    replace_line_in_file(fN, 25, line25, fN)

    #lines 34 and 40
    line34="major_elements<-c("
    line40="bulk_definitions<-c(list(\"{1;1}_{1;1}\"=c("
    for X in Comp:
        strAdd = "\""+X+"\","
        valAdd = "\""+str(Comp[X])+"\","
        line34=line34+strAdd
        line40=line40+valAdd
    line34=line34[:-1]
    line34=line34+")"
    line40=line40+"\"100\")))"
    replace_line_in_file(fN, 34, line34, fN)
    replace_line_in_file(fN, 40, line40, fN)


def convert_rdata_list_to_dict(file_path):
    r = robjects.r
    r.load(file_path)
    r_mat = r['crust'][0][0]
    
    rowNames=r.rownames(r_mat)
    colNames=r.colnames(r_mat)
    
    py_dicts = convert_float_matrix_to_dict(r_mat,colNames, rowNames)
    return py_dicts, rowNames

def convert_float_matrix_to_dict(float_matrix, colnames,rownames):
    ncol=len(colnames)
    nrow=len(rownames)
    data = list(float_matrix)
    #print(rownames)
    #print(data)
    result_dicts = [{colnames[j]: data[j*nrow+i] for j in range(0,ncol)} for i in range(nrow)]
    return result_dicts


def extractRData():
    RDataLoc = fileRData
    dataDict = convert_rdata_list_to_dict(RDataLoc)
    return dataDict


def replace_line_in_file(file_path, line_number, replacement_text, output_file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    if 1 <= line_number <= len(lines):
        lines[line_number - 1] = replacement_text + '\n'

        with open(output_file_path, 'w') as file:
            file.writelines(lines)
    else:
        print(f"Error: Line number {line_number} is out of range.")








