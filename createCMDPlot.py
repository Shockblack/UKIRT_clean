#-----------------------------------------------------------------------
# File name: createCMDPlot.py
# 
# Creates figures of CMD's using the cmd class from ukirtSurot.py. This
# file is meant to handle all the plotting capabilities and needs of CMDs.
# This will take and adapt a lot of the plotting functions from the
# original plotCMD.py file. 
#
# NOTES:
# Requires an external package not default in Anaconda which is an extension
# to matplotlib in order to plot densities efficiently. Package can be found
# at https://github.com/astrofrog/mpl-scatter-density and installed easily
# with "pip install mpl-scatter-density" in the command line.
#
# Programmer: Aiden Zelakiewicz (zelakiewicz.1@osu.edu)
#
# For questions, reach out to:
#   Aiden Zelakiewicz   (zelakiewicz.1@osu.edu)
#   Samson Johnson      (johnson.7080@osu.edu)
#
# Revision History:
#   20-Apr-2022 :   File Created
#   22-Apr-2022 :   Added reddening vector 
#   23-Aug-2022 :   Added argparse for command line arguments
#-----------------------------------------------------------------------

# Gathering all our imports
import matplotlib.pyplot as plt
import createCMD
import parameters as pram
from astropy import coordinates as coord
from astropy import units as u
import numpy as np
import mpl_scatter_density
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
import argparse
from matplotlib import cm
from matplotlib.colors import Normalize 
from scipy.interpolate import interpn

# Create command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--field', type=int, default=39, help='Field index for the cmd')
parser.add_argument('-a', '--altmag', nargs=2, type=float, default=[12,15.5], help='Magnitude bound for fit')
parser.add_argument('-d','--delta', nargs=2, type=float, default=[1.0,2.0], help='Color bound for fit')
opt = parser.parse_args()

def plotDensityCMD(cmd,limit_dict=None, plotvec=False, coeffs=None, savefig=False):
    """Creates a CMD of a given area with the density of stars
    plotted on a color scale. This can be a field of any size, as
    long as it can be passed into the cmd class.

    Parameters
    ----------
    cmd : class
        Class object with all the cmd data of a given area.
    limit_dict : dict, optional
        Dictionary which holds the limits used to determining
        the RC, by default None
    plotvec : bool, optional
        Whether the reddening vector is being calculated.
        Will plot the vector alongside some others if True, by default False
    coeffs : list, optional
        If plotvec=True, the coefficients for the reddening vector so it can be
        plotted. Has the form [slope, intercept]. By default None
    savefig : bool, optional
        Determines whether the figure will be saved or just shown, by default False
    """

    # Creating a custom version of viridis with a white background
    # Credit to user @np8 in post:
    # https://stackoverflow.com/questions/20105364/how-can-i-make-a-scatter-plot-colored-by-density-in-matplotlib
    white_viridis = LinearSegmentedColormap.from_list('white_viridis', [
        (0, '#ffffff'),
        (1e-20, '#440053'),
        (0.2, '#404388'),
        (0.4, '#2a788e'),
        (0.6, '#21a784'),
        (0.8, '#78d151'),
        (1, '#fde624'),
    ], N=256)

    white_rdylbl = LinearSegmentedColormap.from_list('white_rdylbl', [
        (0, '#ffffff'),
        (1e-20, '#4575b4'),
        (0.17, '#91bfdb'),
        (0.34, '#e0f3f8'),
        (0.51, '#ffffbf'),
        (0.68, '#fee090'),
        (0.85, '#fc8d59'),
        (1, '#d73027'),
    ], N=256)

    # Creating the plotting environment
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection='scatter_density')
    # ax = fig.add_subplot(1, 1, 1)
    ax.set_title(r'($l,b$)=(%.2f,%.2f)'%(cmd.l,cmd.b))
    # Creating the density plot using mpl extension package mpl_scatter_density
    density = ax.scatter_density(cmd.filterStarDict['delta'], cmd.filterStarDict['altmag'], cmap = white_rdylbl, dpi=72, vmin=0, vmax=300)
    # density_scatter(cmd.filterStarDict['delta'], cmd.filterStarDict['altmag'], ax=ax, cmap = 'cividis', s=0.1, bins=100)
    fig.colorbar(density, label='Density')

    # Creates a box of the stars used in the reddening vector fit
    if type(limit_dict)!=type(None):
        
        limit_box = []
        keys = []
        for key in limit_dict.keys():
            keys.append(key)
        limit_box.append([limit_dict[keys[0]][0],limit_dict[keys[0]][1],limit_dict[keys[0]][1],limit_dict[keys[0]][0],limit_dict[keys[0]][0]])
        limit_box.append([limit_dict[keys[1]][0],limit_dict[keys[1]][0],limit_dict[keys[1]][1],limit_dict[keys[1]][1],limit_dict[keys[1]][0]])
        
        # plots the box
        ax.plot(limit_box[1],limit_box[0],'C1',linewidth=2)

    # If the reddening vector is calculated, plot it along with Nishiyama et al. 2009 value
    if plotvec == True:
        fitx = np.linspace(cmd.filterStarDict['delta'].min(),cmd.filterStarDict['delta'].max(),100)
        ax.plot(fitx,np.polyval(coeffs,fitx), color='red', label='Derived $R_V={}$'.format(coeffs[0]))
        # ax.plot(cmd.vecFitList['delta'],cmd.vecFitList['altmag'],'.r', label='Fitted Points')
        # ax.plot(fitx,1.46*(fitx-0.4)+13.1, label='Nishiyama et al. 2009')
        ax.plot(fitx,1.5*(fitx-0.4)+13.1, label='Cardelli et al. 1989')

    # Some stylizing
    ax.set_xlabel(r'$\it{H-K}$')
    ax.set_ylabel(r'$\it{K}$')
    ax.legend(loc='upper right')
    
    ax.set_xlim([-.4,2])
    ax.set_ylim([18.5,11])

    # ax.set_xlim([0,2])
    # ax.set_ylim([16,12])

    if savefig==True:
        plt.savefig('../misc_figs/cmd1.pdf')
    else:
        plt.show()


def density_scatter( x , y, ax = None, sort = True, bins = 20, **kwargs )   :
    """
    Scatter plot colored by 2d histogram
    """
    if ax is None :
        fig , ax = plt.subplots()
    data , x_e, y_e = np.histogram2d( x, y, bins = bins, density = True )
    z = interpn( ( 0.5*(x_e[1:] + x_e[:-1]) , 0.5*(y_e[1:]+y_e[:-1]) ) , data , np.vstack([x,y]).T , method = "splinef2d", bounds_error = False)

    #To be sure to plot all data
    z[np.where(np.isnan(z))] = 0.0

    # Sort the points by density, so that the densest points are plotted last
    if sort :
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]

    ax.scatter( x, y, c=z, **kwargs )

    # norm = Normalize(vmin = np.min(z), vmax = np.max(z))
    # cbar = fig.colorbar(cm.ScalarMappable(norm = norm), ax=ax)
    # cbar.ax.set_ylabel('Density')

    return ax

if __name__ == '__main__':
    import csv
    rc_limits_byeye = np.loadtxt('./RC_limits_byEye.csv',delimiter=',', skiprows=1)
    field_data = rc_limits_byeye[opt.field]
    delta = field_data[1:3]
    altmag = field_data[3:5]

    # Creates the cmd_test object and dictionary for vector calculating
    rc_dict={'altmag':opt.altmag, 'delta':opt.delta, 'altMAD':[-0.1,0.1], 'MAD':[-0.1,0.1]}
    # rc_dict={'altmag':altmag, 'delta':delta, 'altMAD':[-0.1,0.1], 'MAD':[-0.1,0.1]}
    for i in range(1):
        i=0
        cmd_test = createCMD.cmd(findvec=True, field_ind=i, fieldType='field',rc_dict=rc_dict)
    
    # Calls the calculation to find the reddening vector
        print("The reddening vector is: " + str(cmd_test.coeffs[0]))
    # Plots the density and vector
        # plotDensityCMD(cmd_test, limit_dict=rc_dict, plotvec=True, coeffs=cmd_test.coeffs, savefig=False)
        plotDensityCMD(cmd_test, limit_dict=None, plotvec=False, coeffs=cmd_test.coeffs, savefig=False)


    # append_bool = input("Do you want to append the vector to the file? (y/n) \n")
    # if append_bool == 'y':
        # with open('RC_limits_byEye.csv', 'a') as csvfile:
            # writer = csv.writer(csvfile)
            # writer.writerow([opt.field, opt.delta[0], opt.delta[1], opt.altmag[0], opt.altmag[1]])

    # with open('red_vecs.csv', 'a') as csvfile:
        # writer = csv.writer(csvfile)
        # writer.writerow([opt.field, cmd_test.coeffs[0]])