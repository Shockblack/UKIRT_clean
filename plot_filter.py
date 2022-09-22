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
#------------------------------------------------------------------------------

# Importing necessary packages and files
import matplotlib.pyplot as plt
import pandas as pd

# Importing the filter data from the SVO project
ukirt = pd.read_csv('../data/filter/UKIRT_WFCAM.K.dat', delimiter=' ', names=['wavelength', 'transmission'])
vvv = pd.read_csv('../data/filter/Paranal_VISTA.Ks.dat', delimiter=' ', names=['wavelength', 'transmission'])
twomass = pd.read_csv('../data/filter/2MASS_2MASS.Ks.dat',delimiter=' ', names=['wavelength','transmission'])

# Modifying the matplotlib settings
plt.style.use('ggplot')

# Plotting the filters
plt.plot(ukirt['wavelength'],ukirt['transmission'],label='UKIRT K')
plt.plot(vvv['wavelength'],vvv['transmission'],label='VVV Ks')
plt.plot(twomass['wavelength'],twomass['transmission'],label='2MASS Ks')

# Labelling the plot
plt.xlabel('Wavelength (Angstroms)')
plt.ylabel('Transmission')
plt.legend()
plt.savefig('figs/filter_comparison.pdf')