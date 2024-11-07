# WOW  ---   Wow! Ocean Worlds.  

## Overview:  
This numerical modeling framework simulates the thermal and compositional evolution of water-rich planetary bodies. It is designed for modeling ocean world interiors, incorporating water-rock reactions, mineral metamorphism, and ice freezing/melting.


## ---- Installation instructions ----

### 1). WOW requires the following software
  
Python with following packages:  
	&emsp; - numpy  
	&emsp; - scipy  
	&emsp; - matplotlib  
	&emsp; - rpy2  
	&emsp; - pandas  
	&emsp; - molmass  
	&emsp; - Jupyter notebook: https://jupyter.org/install  
R: https://www.r-project.org/  
Rcrust: https://www.sun.ac.za/english/faculty/science/earthsciences/rcrust  
Perple_X: https://www.perplex.ethz.ch/  
EQ3/6: https://github.com/LLNL/EQ3_6  

### 2). Clone the WOW repository into a new directory  

### 3). Before the notebook will work, there are file paths that need to be updated.   
Open WOW.py and update the three file paths listed:  
	&emsp; main_directory - This is the location of the WOW repository on your machine.  
	&emsp; rcrust_directory - This is the location of your Rcrust installation on your machine.  
	&emsp; eq36_directory - This is the location of your EQ36 installation on your machine.  

### 4). Move Rcrust and EQ36 project and data files.   
&emsp; - Move the WOW Rcrust project into your Rcrust project folder.  
&emsp; - Move the Rcrust data files from the WOW folder to your Rcrust data file folder.  
&emsp; - Move the EQ36 database files into the db folder of your EQ36 installation.   

### 5). Once all required software has been successfully installed, and all paths updated,  
Open the "Demo" Jupyter notebook  
