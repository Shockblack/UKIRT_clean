#-----------------------------------------------------------------------
# File name: createMap.py
# 
# Pseudo-master file which will call the other files to calculate the
# reddening and extinction for all pixels in all subfields and save them
# to a csv file. 
#
#-----------------------------------------------------------------------

# Adding all necessary imports
import createCMD
import calcReddening
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as pat
from matplotlib import cm
from astropy import coordinates as coord
from astropy import units as u 
import numpy as np
import parameters as pram
import ipdb

def read_map(filename):

    map_data = np.loadtxt(filename,delimiter=',')
    return map_data

def calc_angle(ra1,dec1,ra2,dec2):
    c1 = coord.SkyCoord(ra=ra1*u.degree,dec=dec1*u.degree,frame='icrs')
    c2 = coord.SkyCoord(ra=ra2*u.degree,dec=dec2*u.degree,frame='icrs')
    angle = np.degrees(np.arctan2(c2.galactic.b.degree-c1.galactic.b.degree,c2.galactic.l.degree-c1.galactic.l.degree))
    return angle


def plot_grid_map(map_data, axis=10, cb_label=r'$A(K)$', path='figs/', figname='map.pdf'):

    # Creating the plotting framework
    fig, ax = plt.subplots()
    # Building the color bar
    good_data = map_data[:,axis][~np.isnan(map_data[:,axis])]
    cmap_max = np.max(good_data)
    cmap_min = np.min(good_data)
    cmap_max = np.percentile(good_data,99) 
    cmap_min = np.percentile(good_data,5)
    cmap_min = 0.1
    cmap_max = 2.5

    norm = matplotlib.colors.Normalize(vmin=cmap_min,vmax=cmap_max)
    cmap = cm.get_cmap('cividis_r')
    #sc = ax.scatter(map_data[:,0],map_data[:,1],facecolor=cmap(norm(map_data[:,axis])))
    sc = ax.scatter(map_data[:,2]+1e3,map_data[:,3],marker='.',c=map_data[:,axis],cmap=cmap,vmin=cmap_min,vmax=cmap_max)#facecolor=cmap(norm(map_data[:,axis])))
    cbar = fig.colorbar(sc,ax=ax)#,vmax=cmap_max,vmin=cmap_min)
    cbar.set_label(cb_label)

    pixel = map_data[0]
    angle = calc_angle(pixel[0]-pram.arcmin/60./2.,pixel[1]-pram.arcmin/60./2.,pixel[0]-pram.arcmin/60./2.,pixel[1]+pram.arcmin/60./2.)

    for i in range(len(map_data)):
        pixel = map_data[i]

        if pixel[2]>160:
            rect = pat.Rectangle((pixel[2]-360,pixel[3]),pram.arcmin/60.,pram.arcmin/60.,angle=angle,facecolor=cmap(norm(pixel[axis])))
        else:
            rect = pat.Rectangle((pixel[2],pixel[3]),pram.arcmin/60.,pram.arcmin/60.,angle=angle,facecolor=cmap(norm(pixel[axis])))

        ax.add_patch(rect)

    ax.set_xlim(2.8,-2.2)
    ax.set_ylim(-2.8,2.1)

    ax.set_xlabel('$l$')
    ax.set_ylabel('$b$')

    ax.set_aspect('equal')
    


if __name__=='__main__':
    meta_map_data = []
    map_data = read_map('maps/mcmc_small.map')
    # map_data = read_map('maps/mcmc_map_1.5_Bconst.map')
    
    # map_data[:,6]=map_data[:,6]-12.93
    map_data[:,8]=map_data[:,8]-12.93


    plot_grid_map(map_data,axis=8, cb_label=r'$A(K)$')
    plt.show()