#------------------------------------------------------------------------------
# filename: plot_filter.py
#
# Creates a figure which compares the filters for the UKIRT, VVV, and 2MASS
# surveys in the K and Ks bands.
#
# Acknowledgements:
#   This research has made use of the Spanish Virtual Observatory (https://svo.cab.inta-csic.es)
#   project funded by MCIN/AEI/10.13039/501100011033/ through grant PID2020-112949GB-I00
#
# Programmer: Aiden Zelakiewicz
#
# Revision History:
#   27-Aug-2022 : File Created
#   17-Aug-2023 : Scaling the data to 1
#   21-Aug-2023 : Updated colors
#------------------------------------------------------------------------------

# Importing necessary packages and files
import matplotlib.pyplot as plt
import pandas as pd

# Importing the filter data from the SVO project
ukirt = pd.read_csv('../data/filter/UKIRT_WFCAM.K.dat', delimiter=' ', names=['wavelength', 'transmission'])
vvv = pd.read_csv('../data/filter/Paranal_VISTA.Ks.dat', delimiter=' ', names=['wavelength', 'transmission'])
twomass = pd.read_csv('../data/filter/2MASS_2MASS.Ks.dat',delimiter=' ', names=['wavelength','transmission'])

# Scale the transmission data to 1
ukirt['transmission'] = ukirt['transmission']/max(ukirt['transmission'])
vvv['transmission'] = vvv['transmission']/max(vvv['transmission'])
twomass['transmission'] = twomass['transmission']/max(twomass['transmission'])

# Modifying the matplotlib settings
plt.style.use('az-paper-twocol')
color_list = ['#F27405','#BF0404','k']
# Plotting the filters
plt.plot(ukirt['wavelength'],ukirt['transmission'],label='UKIRT K', c=color_list[0])
plt.plot(vvv['wavelength'],vvv['transmission'],label='VVV Ks', c=color_list[1])
plt.plot(twomass['wavelength'],twomass['transmission'],label='2MASS Ks', c=color_list[2])

# Labelling the plot
plt.xlabel('Wavelength (Angstroms)', fontsize=12)
plt.ylabel('Transmission', fontsize=12)
plt.legend(fontsize=10)
plt.savefig('paperfigs/filters.pdf')
# plt.show()