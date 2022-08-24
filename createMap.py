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
    cmap_max = 1.8#np.percentile(good_data,95) 
    cmap_min = 0.1#np.percentile(good_data,5)

    norm = matplotlib.colors.Normalize(vmin=cmap_min,vmax=cmap_max)
    cmap = cm.get_cmap('cividis_r')
    #sc = ax.scatter(map_data[:,0],map_data[:,1],facecolor=cmap(norm(map_data[:,axis])))
    sc = ax.scatter(map_data[:,2]+1e3,map_data[:,3],marker='.',c=map_data[:,axis],cmap=cmap,vmin=cmap_min,vmax=cmap_max)#facecolor=cmap(norm(map_data[:,axis])))
    cbar = fig.colorbar(sc,ax=ax)#,vmax=cmap_max,vmin=cmap_min)
    cbar.set_label(cb_label)

    pixel = map_data[0]
    angle = calc_angle(pixel[0]-pixel[4]/2.,pixel[1]-pixel[4]/2.,pixel[0]-pixel[4]/2.,pixel[1]+pixel[4]/2.)

    for i in range(len(map_data)):
        pixel = map_data[i]

        if pixel[2]>160:
            rect = pat.Rectangle((pixel[2]-360,pixel[3]),pixel[4],pixel[4],angle=angle,facecolor=cmap(norm(pixel[axis])))
        else:
            rect = pat.Rectangle((pixel[2],pixel[3]),pixel[4],pixel[4],angle=angle,facecolor=cmap(norm(pixel[axis])))

        ax.add_patch(rect)

    ax.set_xlim(2.8,-2.2)
    ax.set_ylim(-2.8,2.1)

    ax.set_xlabel('$l$')
    ax.set_ylabel('$b$')

    ax.set_aspect('equal')
    


if __name__=='__main__':
    meta_map_data = []
    for i in range(56):
        index = str(i)
        map_data = read_map('maps/field_pixel_data_'+index+'.map')
        for star in map_data:
            meta_map_data.append(star)

    meta_map_data = np.array(meta_map_data)
    ipdb.set_trace()
    plot_grid_map(meta_map_data,axis=10)
    plt.show()