setwd("D:/Rcrust/code")
  manual_load<<-function(working_file,projects_directory=paste0(substring(getwd(),1,nchar(getwd())-4),"Projects")){
    source(paste0(projects_directory,"/",working_file,"/Inputs/",working_file,".txt"))
    source("main.r")
  }