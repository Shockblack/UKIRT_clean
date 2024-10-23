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
plt.style.use('az-paper-twocol')
# plt.rcParams.update({
#     "text.usetex": True,
#     "font.family": "Computer Modern Roman"
# })

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
    fig, ax = plt.subplots(figsize=(8,6))
    ax.fill_between([-3,3],-0.8,0.8,alpha=0.1,color='k')
    # Building the color bar
    good_data = map_data[:,axis][~np.isnan(map_data[:,axis])]
    cmap_max = np.max(good_data)
    cmap_min = np.min(good_data)
    cmap_max = np.percentile(good_data,99) 
    cmap_min = np.percentile(good_data,1)
    cmap_min = 0.1
    cmap_max = 2.5

    norm = matplotlib.colors.Normalize(vmin=cmap_min,vmax=cmap_max)
    cmap = cm.get_cmap('gist_heat_r')
    #sc = ax.scatter(map_data[:,0],map_data[:,1],facecolor=cmap(norm(map_data[:,axis])))
    sc = ax.scatter(map_data[:,2]+1e3,map_data[:,3],marker='.',c=map_data[:,axis],cmap=cmap,vmin=cmap_min,vmax=cmap_max)#facecolor=cmap(norm(map_data[:,axis])))
    cbar = fig.colorbar(sc,ax=ax, fraction=0.046, pad=0.04)#,vmax=cmap_max,vmin=cmap_min)
    cbar.set_label(cb_label, fontsize=16)
    cbar.ax.tick_params(labelsize=16)

    pixel = map_data[0]
    angle = calc_angle(pixel[0]-pram.arcmin/60./2.,pixel[1]-pram.arcmin/60./2.,pixel[0]-pram.arcmin/60./2.,pixel[1]+pram.arcmin/60./2.)

    for i in range(len(map_data)):
        pixel = map_data[i]
        l = pixel[2]
        b = pixel[3]

        if l > 180:
            l = l-360

        rect = pat.Rectangle((l,b),pram.arcmin/60.,pram.arcmin/60.,angle=angle,facecolor=cmap(norm(pixel[axis])))

        ax.add_patch(rect)
        # if l > -0.2 and l < 0.08 and b > -0.08 and b <0.088:
            # ax.text(l,b,str(i),color='k',fontsize=8)

    ax.set_xlim(2.8,-2.2)
    ax.set_ylim(-2.8,2.1)

    ax.set_xlabel('Galactic Longitude ($l$)', fontsize=16)
    ax.set_ylabel('Galactic Latitude ($b$)', fontsize=16)
    ax.hlines(0.8,-3,3, alpha=0.5, color='k',linestyle='--')
    ax.hlines(-0.8,-3,3, alpha=0.5, color='k',linestyle='--')
    

    ax.set_aspect('equal')
    plt.tight_layout()
    ax.tick_params(axis='both',which='both',direction='in',top=True,right=True, labelsize=16)
    


if __name__=='__main__':
    meta_map_data = []
    map_data = read_map('paperdata/mcmc_map.map')
    # map_data = read_map('maps/nm_map_2_B.map')
    # map_data = read_map('maps/mcmc_map_1.5_Bconst.map')
    
    # map_data[:,6]=map_data[:,6]-12.94
    axis = 8
    map_data[:,axis]=map_data[:,axis] - 12.97 # 12.94 for K, 12.93 for Ks
    # map_data[:,axis]=map_data[:,axis] - 0.16
    # axis = 5
    # mad = np.nanmedian(abs(map_data[:,5]-np.nanmedian(map_data[:,5])))
    
    # map_data[:,5]=abs(map_data[:,5]-np.nanmedian(map_data[:,5]))/mad
    plot_grid_map(map_data,axis=axis, cb_label=r'$A(K)$')
    # fields = []
    # layout = 'fieldLocations/romanfields.txt'
    # from itertools import groupby
    # with open(layout) as fp:
    #     for k, g in groupby(fp, lambda x: x.startswith(' ')):
    #         if not k:
    #             fields.append(np.array([[float(x) for x in d.split()] for d in g if len(d.strip())]))

    # for f in fields:
    #     plt.plot(f[:,1],f[:,2],'k-',lw=3)
    # plt.savefig('paperfigs/map.pdf')
    plt.show()