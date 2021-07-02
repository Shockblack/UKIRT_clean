#Uses input in order to streamline map production
#Trying to make this file a "one stop shop" for using the make_map.py file
#edit: Now including mapper.py (WOW, such value!)
from sys import exit
import make_map as mm
import mapper as mpp
import ipdb
import parameters as pram

print(
"""
Entering map creator. Please enter parameters when prompted.
If no parameter is specified, error will be raised (default feature tbd).
WARNING: Make sure you created the .map file for the desired year and edge length within a /maps folder.
"""
    )

print(
"""
Please input selected year:
====================================
| 2015 | 2016 | 2017 | 2018 | 2019 |
====================================
"""
)

yearIn = int(float(input('Enter Year: ')))
if yearIn == 2015 or yearIn == 2016 or yearIn == 2017 or yearIn == 2018 or yearIn == 2019:
    print("Input accepted! Selected year:",yearIn)
    print("")
else:
    print(
"""
===========================
ERROR: Input not valid year
===========================
"""
    )
    exit()
try:
    edgeIn = float(input('Please input the desired edge length in arcminutes: '))
    print("Input accepted! Selected edge length of",edgeIn,"arcminutes")
except:
    print(
"""
===================================
ERROR: Input not a float or integer
===================================
"""
    )
    exit()

print(
"""
Do you need to create a .map file? Please enter
yes (y) or no (n).
"""
) 
mapIn = input("[Yes/No]: ")
mapIn = mapIn.lower()[:1]
if mapIn == 'y':
    test_map = mpp.map(year=yearIn,edge_length=float(edgeIn)/60)#,ra_lims=[269,280],dec_lims=[-31,-26])#,ra_lims=[268.75,289.25],dec_lims=[-29.75,-29.25])
    test_map.get_fields_stats()
    test_map.gen_grid()
    test_map.filter_grid()

    mapnameIn = 'map_'+pram.phot+'_'+str(yearIn)+'_'+str(edgeIn)
    test_map.fit_map(filebase=mapnameIn)
    test_map.save_map(filename=mapnameIn)
elif mapIn == 'n':
    print("""Using existing .map file...
    """)
else:
    print(
"""
========================
ERROR: Not a valid input
========================
"""
    )
    exit()


print(
"""Input number corresponding to map type:
=========================================================
| 1. N_stars     | 2. A     | 3. Aerr     | 4. B        |
| 5. Berr        | 6. M_RC  | 7. M_RCerr  | 8. sigma_RC |
| 9. sigma_RCerr | 10. N_RC | 11. N_RCerr |             |
========================================================="""
)
#ipdb.set_trace()
axisHeader = ['RA','DEC','l','b','edgelength','N_stars','A','Aerr','B','Berr','M_RC','M_RCerr',\
    'sigma_RC','sigma_RCerr','N_RC','N_RCerr','FinalColorRC','RealSTDcolorOpt','STDtotal','VarianceMin','MU1opt','MU2opt']

try:
    axisIn = int(float(input('Enter Integer 1-11: ')))
    print("Input accepted! Selected type:",axisIn)
    #Labels for the color bars
    label=[r'$N-Stars$',r'$A$',r'$A-error$',r'$B$',r'$B-error$',r'$A(k)$',r'$A(k)-error$',u'$sigma_{RC} $',r'$sigma_{RC}-error $',r'$N_{RC}$',r'$N_{RC}-error$',r'FinalColorRC',r'RealSTDcolorOpt']
    cblabel=label[axisIn-1]
    #Aligning the value inputed with the location in the list
    axisIn += 4
except:
    print(
"""
===================================
ERROR: Input not a float or integer
===================================
"""
    )
    exit()

#Basically just determining whether 13.1 needs to be subtracted
#from M_RC in order to make a reddening/A(K) figure
if axisIn == 10:
    funcIn = mm.A_K
#elif axisIn == 14:
#    funcIn = mm.N_RC
else:
    funcIn = mm.do_nothing


print(
"""
Please input the selected file name for the outputted map.
Make sure to exclude any file tags as .pdf is assumed. Also make
sure a 'figs' folder is included within the project as that is the
defaulted output path. It default contains the string below. Any
input will add it onto the end.
mtype_photom_year_arcmin
"""
)
nameadd = input('Filename: ')
if axisIn == 10:
    mtype = "AK"
elif axisIn == 11:
    mtype = "AKerr"
else:
    mtype = axisHeader[axisIn]

if bool(nameadd) == True:
    fignameIn = mtype+'_'+pram.phot+'_'+str(yearIn)+'_'+str(int(edgeIn))+'_'+str(nameadd)
else:
    fignameIn = mtype+'_'+pram.phot+'_'+str(yearIn)+'_'+str(int(edgeIn))
print(
"""
Creating Map...
"""
)

mapname = mm.read_map('maps/map_'+pram.phot+'_'+str(yearIn)+'_'+str(int(edgeIn))+'.map')
mm.plot_grid_map(mapname,func=funcIn,cb_label=cblabel,axis=axisIn,figname=fignameIn+'.pdf')

print(
"""
Figure successfully created and saved to figs directory!
"""
)