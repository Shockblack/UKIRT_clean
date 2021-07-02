#File for inputing global setting changes to be used within the meat of the code.
#Allows for streamlining of changes across multiple files

#defining the arcminutes used for the run (pixel size is arcmin/60)
arcmin = 2

year = 2017
#year = 2018
#year = 2019

phot = 'PSF'
#phot = 'CASU'

if year == 2019: #I dont remember what this was for...
    pass
elif year == 2017:
    pass
#Method for IC's... pix might be removed eventually
#method = 'pix'
method = 'mag'

#Default pathways
