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

#Default pathways... Please enter the path to the project folder here so no internal changes need to be made!
projpath = '/mnt/a/documents/files/surp/ukirt-exinction-maps'

brown = '#260101'
red = '#BF0404'
orange = '#F27405'

map_header = ['ra', 'dec', 'l', 'b', 'A', 'A_err_low', 'A_err_up', 'B', 'B_err_low', 'B_err_up', \
                'MRC', 'MRC_err_low', 'MRC_err_up', \
                'SIGMA', 'SIGMA_err_low', 'SIGMA_err_up', \
                'NRC', 'NRC_err_low', 'NRC_err_up', \
                'EWRC', 'EWRC_err_low', 'EWRC_err_up', \
                'HK', 'sigmaHK'
            ]

ext_header = ['l', 'b', 'AK', 'AK_err_low', 'AK_err_up', 'EHK', 'EHK_err']