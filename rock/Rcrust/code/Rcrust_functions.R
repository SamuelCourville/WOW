#Repository for Rcrust functions
crust_to_gcdkit<-function(crust,choose_columns=NULL,choose_rows=NULL,choose_points="All",GCDkitGUI=FALSE,source.first=TRUE,source.plugins=TRUE){
#crust_to_gcdkit(crust,choose_points="{1;1}")
#sends crust object to gcdkit
#note this function has to assign globally so that GCDkit can see it
	data_crust<-as.data.frame(data_file(crust,x_n,y_n,choose_columns,choose_rows,choose_points))
	rownames(data_crust)<-data_crust[,1]
	require(GCDkitDevelop)
	#library(GCDkit)
		stopApp(returnValue = invisible())
	#source(paste(gcdx.dir,"/GCDkit.r",sep=""))
	#setwd(gcdx.dir)
	accessVar(data_crust,source.first=source.first,source.plugins=source.plugins,poke.data=FALSE,GUI=GCDkitGUI)
	print("To launch the Rcrust GUI type Rcrust() into the console and press enter")
	  #Provide again the launch with GUI function in order to get back
  Rcrust<<-function(){
   #If working directory is x\Projects\y then set to x\code
    if(length(grep("Rcrust/Projects/",getwd()))==1){
	setwd(paste0(strsplit(getwd(),split="Projects")[[1]][1],"code"))
	}
    runApp()
  }
}

assign_label<-function(crust,from_label,to_label,label_name,label_value){
#adds column to crust and populates it with value over range
#example
# label_name<-"Protolith"
# from_label<-"{1;1}"
# to_label<-"{1;1}"
# label_value<-"Sample1"
a<-unlist(strsplit(gsub("\\{","",gsub("\\}","",from_label)),split=";"))
b<-unlist(strsplit(gsub("\\{","",gsub("\\}","",to_label)),split=";"))
for(x_i in a[1]:b[1]){
for(y_i in a[2]:b[2]){
if(!any(colnames(crust[[y_i]][[x_i]])==label_name)){
lab_column<-matrix(label_value,nrow(crust[[y_i]][[x_i]]),1)
colnames(lab_column)<-label_name
crust[[y_i]][[x_i]]<-cbind(crust[[y_i]][[x_i]],lab_column)
}else{
crust[[y_i]][[x_i]][,lab_column]<-label_value
}
}
}
return(crust)
}

data_file<-function(crust,x_n=length(crust[[1]]),y_n=length(crust),choose_columns=NULL,choose_rows=NULL,choose_points="All"){
  #######################################################################
  #Outputs select_data list
  #Settings for outputing data_file

  # choose_columns    options =   from crust ("wt%",comps,"mass","V(J/bar)","H(J)","Gruneisen_T","Ks(bar)","Mu(bar)","V0(km/s)","Vp(km/s)","Vs(km/s)","Vp/Vs","Rho(kg/m3)","Cp(J/K)","alpha(1/K)","beta(1/bar)","S(J/K)","N(g)","Cp/Cv")
  #                               from other ("Phase","y_i","x_i","Pressure(kbar)","Temperature(C)")
  #                               Brief
  #                               All
  #                   default = All
  #choose_columns<-c("All")
  #choose_columns<-c("Brief")

  #choose_rows       default = All
  #choose_rows       options =  Reactive subsystem
  #                             Extract subsystem
  #                             Addition subsystem
  #choose_rows<-c("All")
  #choose_rows<-c("Bio(TCC)_rs","Bulk_rs")

  #choose_points           default = All

  #######################################################################
  #fix-tag: number duplicate names (feldspar, mica)
  validate(need(!choose_points=="","To select points enter arguments seperated by commas of the form {x_a;y_a} for single points or {x_a;y_a}:{x_b;y_b} for ranges where a<=b<=n"))
  if(choose_points=="All"){choose_points<-paste("{1;1}:{",x_n,";",y_n,"}",sep="")}
  choose_points<-unlist(strsplit(choose_points,split=","))
  choose_points<-unlist(strsplit(gsub("\\{","",gsub("\\}","",choose_points)),split=","))
  choose_points<-strsplit(choose_points,split=":|;")
  data_out<-NULL
  for(i in 1:length(choose_points)){
  x_a<-as.numeric(choose_points[[i]][1])
  y_a<-as.numeric(choose_points[[i]][2])
  if(length(choose_points[[i]])>2){
  x_b<-as.numeric(choose_points[[i]][3])
  y_b<-as.numeric(choose_points[[i]][4])
  }else{
  x_b<-as.numeric(choose_points[[i]][1])
  y_b<-as.numeric(choose_points[[i]][2])
  }
  #error validation on selection (range must be possible i.e. b>=a,b<=n)
  validate(need(all(x_a<=x_b,y_a<=y_b,x_b<=x_n,y_b<=y_n),"To select points enter arguments seperated by commas of the form {x_a;y_a} for single points or {x_a;y_a}:{x_b;y_b} for ranges where a<=b<=n"))
  for(y_i in y_a:y_b){
    for(x_i in x_a:x_b){
	#fix-tag: create validation here for erronous pixels
	#error validation for point existence
	if(is.null(crust[[y_i]][[x_i]])){
	#fix-tag: only works if point 1;1 is populated - dont know how many columns otherwise
	new_pnt<-matrix(NA,1,(ncol(crust[[1]][[1]])+6))
	new_pnt[1,1]<-paste("Blank",y_i,x_i,sep="_")
	}else{
	  ID<-matrix(paste(rownames(crust[[y_i]][[x_i]]),y_i,x_i,sep="_"),ncol=1)
      colnames(ID)<-"ID"
      phase<-matrix(rownames(crust[[y_i]][[x_i]]),ncol=1)
      colnames(phase)<-"Phase"
      pnt<-matrix(c(y_i,x_i),nrow=nrow(crust[[y_i]][[x_i]]),ncol=2,byrow=TRUE)
      colnames(pnt)<-c("y_i","x_i")
      if(exists("input_pt")){
      p_t<-matrix(input_pt[[y_i]][[x_i]],nrow=nrow(crust[[y_i]][[x_i]]),ncol=2,byrow=TRUE)
      }else{
      p_t<-matrix(0,nrow=nrow(crust[[y_i]][[x_i]]),ncol=2,byrow=TRUE)
      }
      colnames(p_t)<-c("Pressure(kbar)","Temperature(C)")
	  #new_pnt<-cbind(ID,phase,pnt,p_t,signif(crust[[y_i]][[x_i]],4))
	  new_pnt<-cbind(ID,phase,pnt,p_t,crust[[y_i]][[x_i]])
	  rownames(new_pnt)<-ID
	  }
	  data_out<-rbind(data_out,new_pnt)
	}
	}
  }
   if(is.null(choose_columns)){
	choose_columns<-colnames(data_out)
  }
  if(any(choose_columns=="Brief")){
    choose_columns<-union(choose_columns[-which(choose_columns=="Brief")],c("ID","Phase","y_i","x_i","Pressure(kbar)","Temperature(C)","wt%",comps,"mass"))
  }
    if(is.null(choose_rows)){
      select_rows<-1:nrow(data_out)
    }else{
	if(choose_rows=="All"){
	select_rows<-1:nrow(data_out)}else{
      select_rows<-NULL
      if(!is.na(match("Reactive subsystem",choose_rows))){
      rs_rows<-grep("_rs",rownames(data_out))
      choose_rows<-choose_rows[-match("Reactive subsystem",choose_rows)]
      select_rows<-union(select_rows,rs_rows)
      }
      if(!is.na(match("Extract subsystem",choose_rows))){
        rs_rows<-grep("_es",rownames(data_out))
        choose_rows<-choose_rows[-match("Extract subsystem",choose_rows)]
        select_rows<-union(select_rows,rs_rows)
      }
      for(row_arg in choose_rows){
	  chk_names<-NULL
	  for(i in 1:length(rownames(data_out))){
	  chk_names<-c(chk_names,paste(strsplit(rownames(data_out),"_")[[i]][c(-3,-4)],collapse="_"))
	  }
		select_rows<-union(select_rows,which(chk_names==row_arg))
      }
	  }
    }
	data_out<-data_out[sort(select_rows),choose_columns,drop=FALSE]
	rownames(data_out)<-NULL
  return(data_out)
}
	