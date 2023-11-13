# Test script for GCDmodel
setwd("C:/Users/Matthew/Documents/Work/2019-01 Trace element routine")
library(GCDkit)
library(GCDmodel)
kd.in<-read.table(file="yak.kd",sep="\t")
src<-read.table(file="yakC0.txt",sep="\t")
nn<-colnames(src)
src<-as.numeric(src)
names(src)<-nn
liqs<-read.table(file="test_liq",sep="\t")
liqdf<-liqs[1,] # a data.frame
nn<-colnames(liqdf)
liqv<-as.numeric(liqdf) # a vector
names(liqv)<-nn
###
# Major elements
mjr<-c("SiO2","TiO2","Al2O3","FeO","MgO","CaO","Na2O","K2O")
## Cleaning the Kd file
p<-colnames(liqs)[14:21]
#p<-c("Kfs","Bio","Gt_1","Gt_2","Osm","Zrn")
kd.ppx<-ppxCleanKdTable(kd.in,ppxPhases=p)
## batch
bpm<-BatchPM(kd=kd.ppx,
        c0=src,
        pm=liqdf[,"FF"],
        min.props=liqdf[,p],
        cmins=matrix(),melt.arg=list(),dont=character(0))

bpm<-BatchPM(kd=kd.ppx,
             c0=src,
             pm=liqv["FF"],
             min.props=liqv[p],
             cmins=matrix(),melt.arg=list(),dont=character(0))


## Zrn sat
mz<-correctZrnMnzSat(kd=kd.ppx,
                  c0=src,
                  pm=liqv["FF"],
                  min.props=liqdf[,p],
                  melt.arg=list(TT=liqdf[,"Temp"]+273.15,
                                mjrs=liqdf[,mjr],
                                trc=bpm$cL,
                                H2O=liqdf[,"H2O"]
                  ),
                  cmins=bpm$cmins,dont=character(0))


zz<-correctZrnSat(kd=kd.ppx,
              c0=src,
              pm=liqv["FF"],
              min.props=liqdf[,p],
              melt.arg=list(TT=liqdf[,"Temp"]+273.15,
                            mjrs=liqdf[,mjr],
                            trc=bpm$cL
                            ),
              cmins=bpm$cmins,dont=character(0))

mm<-correctMnzSat(kd=kd.ppx,
                  c0=src,
                  pm=liqv["FF"],
                  min.props=liqdf[,p],
                  melt.arg=list(TT=liqdf[,"Temp"]+273.15,
                                mjrs=liqdf[,mjr],
                                H2O=liqdf[,"H2O"],
                                trc=bpm$cL
                  ),
                  cmins=bpm$cmins,dont=character(0))


srcNO<-src
srcNO["Zr"]<-20
bpmNO<-BatchPM(kd=kd.ppx,
             c0=srcNO,
             pm=liqv["FF"],
             min.props=liqv[p],
             cmins=matrix(),melt.arg=list(),dont=character(0))
zzNO<-correctZrnSat(kd=kd.ppx,
                  c0=srcNO,
                  pm=liqv["FF"],
                  min.props=liqdf[,p],
                  melt.arg=list(TT=liqdf[,"Temp"],
                                mjrs=liqdf[,mjr],
                                trc=bpmNO$cL
                  ),
                  cmins=matrix(),dont=character(0))

				  
				  
				  
				  
#calc phases to liqdf (only need liqdf and src to calculate)
liqdf[,mjr]<-calc_phases["Melt",toupper(mjr)]
liqdf[,"FF"]<-calc_phases["Melt","wt%"]
data.frame(liqdf[c(-14:-21)],t(calc_phases[c(-1,-6),"wt%"]))
#source is c0, need to add trace elements to initial comps allowed
?BatchPM


#throw warming if phase proportion is present


GCDkit interfacing
accessVar("test")
#look into new=false etc