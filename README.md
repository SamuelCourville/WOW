# WOW

Wow! Ocean Worlds.  

Overview:  
This numerical modeling framework simulates the thermal and compositional evolution of water-rich planetary bodies. It is designed for modeling ocean world interiors, incorporating water-rock reactions, mineral metamorphism, and ice freezing/melting.


Installation instructions:

WOW requires the following software:  
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

Clone the WOW repository into a new directory  

Before the notebook will work, there are paths that need to be updated.   
Update the commented paths in the following files:  
	- lookupProps.py  
	- Perple_X path  
	- plotHelper.py  
	- eq36.py  

move the Rcrust WOW project into your Rcrust project folder.   

Once all required software have been successfully installed, and all paths updated,  
use the demonstration Jupiter notebook  
