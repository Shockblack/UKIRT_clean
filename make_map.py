import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as pat
from matplotlib import cm
from collections import defaultdict
import ipdb
from readData import readUKIRTfields
from readData import readAltbandCrossCheck
import astropy
from astropy import coordinates as coord
from astropy import units as u
import plotly.express as px

def read_map(filename):
    #ipdb.set_trace()
    map_data = np.loadtxt(filename,delimiter=',')
    return map_data

def plot_point_map(map_data,axis=6):
    fig, ax = plt.subplots()
    
    #build color bar

    good_data = map_data[:,axis][~np.isnan(map_data[:,axis])]
    cmap_max = np.max(good_data)
    cmap_min = np.min(good_data)
    cmap_max = np.percentile(good_data,95)
    cmap_min = np.percentile(good_data,5)

    
    norm = matplotlib.colors.Normalize(vmin=cmap_min,vmax=cmap_max)
    cmap = matplotlib.cm.get_cmap('cividis')

    sc = ax.scatter(map_data[:,0],map_data[:,1],marker='.',c=map_data[:,axis],cmap=cmap,vmin=cmap_min,vmax=cmap_max)#facecolor=cmap(norm(map_data[:,axis])))
    fig.colorbar(sc,ax=ax)#,vmax=cmap_max,vmin=cmap_min)

    ax.set_aspect('equal')
    #plt.show()


def calc_angle(ra1,dec1,ra2,dec2):
    c1 = coord.SkyCoord(ra=ra1*u.degree,dec=dec1*u.degree,frame='icrs')
    c2 = coord.SkyCoord(ra=ra2*u.degree,dec=dec2*u.degree,frame='icrs')
    angle = np.degrees(np.arctan2(c2.galactic.b.degree-c1.galactic.b.degree,c2.galactic.l.degree-c1.galactic.l.degree))
    return angle

def do_nothing(array):
    return array 

def A_K(array):
    M_RC = 13.1
    return array - M_RC

def N_RC(array):
    return


def plot_grid_map(map_data,lb=True,func=do_nothing,axis=7,cb_label=r'$A(K)$',path='figs/',figname='map.pdf'):
    fig, ax = plt.subplots()
    
    #build color bar
    good_data = func(map_data[:,axis][~np.isnan(map_data[:,axis])])
    cmap_max = np.max(good_data)
    cmap_min = np.min(good_data)
    cmap_max = np.percentile(good_data,95) 
    cmap_min = np.percentile(good_data,5)

    norm = matplotlib.colors.Normalize(vmin=cmap_min,vmax=cmap_max)
    cmap = cm.get_cmap('cividis_r')
    #sc = ax.scatter(map_data[:,0],map_data[:,1],facecolor=cmap(norm(map_data[:,axis])))
    sc = ax.scatter(map_data[:,2]+1e3,map_data[:,3],marker='.',c=func(map_data[:,axis]),cmap=cmap,vmin=cmap_min,vmax=cmap_max)#facecolor=cmap(norm(map_data[:,axis])))
    cbar = fig.colorbar(sc,ax=ax)#,vmax=cmap_max,vmin=cmap_min)
    cbar.set_label(cb_label)
    
    #Check location of grids where line discontinuity occurs
    #We first find location then try to find avg mag difference
    

        
    if lb:
        #prep for l,b
        pixel = map_data[0]
        angle = calc_angle(pixel[0]-pixel[4]/2.,pixel[1]-pixel[4]/2.,pixel[0]-pixel[4]/2.,pixel[1]+pixel[4]/2.)
        

        if pixel[2] < 358 and pixel[2]>160:
            print(pixel[0],pixel[1],pixel[2],pixel[3])

        for pixel in map_data:
            #if pixel[3] < -2:
            #    print(pixel[0],pixel[1],pixel[2],pixel[3])
            if np.abs(pixel[3]) <1 and pixel[axis]<13.5:
                print(pixel[0],',',pixel[1],pixel[2],pixel[3])

            if pixel[2]>160:
                rect = pat.Rectangle((pixel[2]-360,pixel[3]),pixel[4],pixel[4],angle=angle,facecolor=cmap(norm(func(pixel[axis]))))
            else:
                rect = pat.Rectangle((pixel[2],pixel[3]),pixel[4],pixel[4],angle=angle,facecolor=cmap(norm(func(pixel[axis]))))
            
            #Lines below for coloring red spots on sharp line changes (UKIRT stuff)
            #if pixel[2]>358.96 and pixel[2] < 359.0113 and pixel[3] > -1.52 and pixel[3] < -1.435:
            #    rect = pat.Rectangle((pixel[2]-360,pixel[3]),pixel[4],pixel[4],angle=angle,facecolor='red')
                
            ax.add_patch(rect)


        ax.set_xlim(2.8,-2.2)
        ax.set_ylim(-2.8,2.1)
        ax.set_xlabel('$l$')
        ax.set_ylabel('$b$')
    
    ax.set_aspect('equal')
    #plt.savefig(path+figname)
    plt.show()



if __name__=='__main__':

    test_map = read_map('maps/map_PSF_2017_2.map')
    
    while True:
        ipdb.set_trace()
        plot_grid_map(test_map,func=A_K,axis=10,figname='MagForUKIRT_PSF_2017.pdf')#[:1000])
        #plot_grid_map(test_map,func=A_K,axis=7)#[:1000])
