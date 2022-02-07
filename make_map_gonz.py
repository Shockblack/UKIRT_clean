import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as pat
from matplotlib import cm
#from matplotlib import cm
#from collections import defaultdict
import ipdb
from readData import readUKIRTfields
from readData import readAltbandCrossCheck
import astropy
from astropy import coordinates as coord
from astropy import units as u
from itertools import groupby



def calc_angle(ra1,dec1,ra2,dec2):
    c1 = coord.SkyCoord(ra=ra1*u.degree,dec=dec1*u.degree,frame='icrs')
    c2 = coord.SkyCoord(ra=ra2*u.degree,dec=dec2*u.degree,frame='icrs')
    
    angle = np.degrees(np.arctan2(c2.galactic.b.degree-c1.galactic.b.degree,c2.galactic.l.degree-c1.galactic.l.degree))
    #angle = np.degrees(np.arctan2(c2.galactic.b.radians-c1.galactic.b.radians,c2.galactic.l.radians-c1.galactic.l.radains))
    #angle = np.degrees(np.arctan2(np.radians(c2.galactic.l.degree-c1.galactic.l.degree),np.radians(c2.galactic.b.degree-c1.galactic.b.degree)))

    return angle

def do_nothing(array):
    return array 

def A_K(array):
    M_RC = 13.1
    return array - M_RC

def E_H_K(array):
    RCcolor = 0.09
    return array - RCcolor

def N_RC(sigma_RC,N_RC):
    sqrt2pi = np.sqrt(2*np.pi)
    vals = np.sqrt(np.pi/(1/2/sigma_RC**2))*N_RC/(sqrt2pi*sigma_RC)
    return vals

class mapper:
    def __init__(self,filename):
        self.read_map(filename)
        try:
            #ipdb.set_trace()
            self.read_gonz(filename.split('.')[0]+'.gonz')
        except:
            print('no gonz data')
        #get angle offset
        pixel = self.map_data[0]
        self.angle = calc_angle(pixel[0]-pixel[4]/2.,pixel[1]-pixel[4]/2.,pixel[0]-pixel[4]/2.,pixel[1]+pixel[4]/2.)
        print(self.angle)
        return None

    def read_gonz(self,filename,path='maps/gonz_maps/'):
        self.gonz_data = np.loadtxt(path+filename,delimiter=' ')
        
    def plot_R(self):
        val1 = A_K(self.map_data[:,10]) 
        val2 = E_H_K(self.map_data[:,16]) 
        self.gen_plot_grid_map(val1/val2,cb_label=r'$R_{HK}$')

    def map_frac_H_K(self,cb_lims=[.1,.5],cb_label=r'$\sigma_{(H-K)}/E(H-K)$',plotname='frac.pdf',plotoutlines=False):
        tmap.gen_plot_grid_map((tmap.map_data[:,17])/(tmap.map_data[:,16]-.09),cb_lims=cb_lims,cb_label=cb_label,plotname=plotname,plotoutlines=plotoutlines)

    def map_sigma_H_K(self,cb_lims=[.03,.4],cb_label='$\sigma_{(H-K)}$',plotname='diff_red.pdf',plotoutlines=False):
        tmap.gen_plot_grid_map(tmap.map_data[:,17],cb_lims=cb_lims,cb_label=cb_label,plotname=plotname,plotoutlines=plotoutlines)

    def map_H_K(self,cb_lims=[.3,1.7],cb_label='$(H-K)$',plotname='color.pdf',plotoutlines=False):
        tmap.gen_plot_grid_map(tmap.map_data[:,16],cb_lims=cb_lims,cb_label=cb_label,plotname=plotname,plotoutlines=plotoutlines)

    def map_E_H_K(self,cb_lims=[.2,1.6],cb_label='$E(H-K)$',plotname='reddening.pdf',plotoutlines=False):
        tmap.gen_plot_grid_map(tmap.map_data[:,16]-.09,cb_lims=cb_lims,cb_label=cb_label,plotname=plotname,plotoutlines=plotoutlines)

    def map_A_K(self,cb_lims=[.0,2.5],cb_label='$A(K)$',plotname='A_K.pdf',plotoutlines=False):
        tmap.gen_plot_grid_map(tmap.map_data[:,10]-13.1,cb_lims=cb_lims,cb_label=cb_label,plotname=plotname,plotoutlines=plotoutlines)

    def map_R_K(self,cb_lims=[.0,2.5],cb_label='$A(K)/E(H-K)$',plotname='A_K.pdf',plotoutlines=False):
        tmap.gen_plot_grid_map((tmap.map_data[:,10]-13.1)/(tmap.map_data[:,16]-.09),cb_lims=cb_lims,cb_label=cb_label,plotname=plotname,plotoutlines=plotoutlines)

    def plot_N_RC(self):
        val1 = self.map_data[:,14]
        #val1 = N_RC(self.map_data[:,12],self.map_data[:,14]) 
        self.gen_plot_grid_map(val1,cb_label=r'$N_RC$')

    def sigma_H_K_hist(self):
        bins = np.linspace(0,.3,100)
        #bin_centers = (bins[1:]+bins[:-1])/2

        #hist,bins = np.histogram(np.log10(self.map_data[:,11],self.map_data[:,10]),bins)
        
        fig,ax = plt.subplots()
        ax.hist((self.map_data[:,17]),bins)
        ax.set_xlabel(r'$\sigma_{(H-K)}$')
        ax.set_ylabel('Number of Pixels')
        ax.set_xlim(0,.2)
        plt.savefig('int_rep_2/sigma_H_K_hist.pdf')
        plt.show()

    def err_M_RC_hist(self):
        bins = np.linspace(-3,-1,100)
        bin_centers = (bins[1:]+bins[:-1])/2

        #hist,bins = np.histogram(np.log10(self.map_data[:,11],self.map_data[:,10]),bins)
        
        fig,ax = plt.subplots()
        ax.hist(np.log10(self.map_data[:,11]/self.map_data[:,10]),bins)
        ax.set_xlabel(r'$\log\left(\sigma_{K_{RC}}/K_{RC}\right)$')
        ax.set_ylabel('Number of Pixels')
        ax.set_xlim(-2.8,-1.25)
        plt.savefig('int_rep_2/A_K_errhist.pdf')
        plt.show()
        

    def plot_ave_stars(self):
        ipdb.set_trace()
        bins = np.linspace(-2.8,2.8,100)
        bin_centers = (bins[1:]+bins[:-1])/2
        
        n_rc_mask = ~np.isnan(self.map_data[:,14])
        s_rc_mask = ~np.isnan(self.map_data[:,12])

        mask = np.logical_and(n_rc_mask,s_rc_mask)
        
        b = self.map_data[:,3][mask]
        n_rc = self.map_data[:,14][mask]
        e_rc  = self.map_data[:,15][mask]
        s_rc = self.map_data[:,12][mask]

        bmask = np.histogram(b, bins)[0] !=0


        aves = (np.histogram(b, bins, weights=N_RC(s_rc,n_rc))[0] / np.histogram(b, bins)[0])
        
        aves1 = (np.histogram(b, bins, weights=n_rc)[0] / np.histogram(b, bins)[0])
        errs = np.sqrt((np.histogram(b, bins, weights=e_rc**2)[0])) / (np.histogram(b, bins)[0])

        ipdb.set_trace()

        fig, ax = plt.subplots()

        #ax.errorbar(bin_centers,np.log10(aves),np.log10(np.sqrt( np.histogram(self.map_data[:,3], bins)[0])),linestyle='None',marker='o')
        #ax.plot(bin_centers,np.log10(aves),'o',label='old')
        #ax.plot(bin_centers,np.log10(aves1),'o',label='new')

        #ax.plot(bin_centers,np.log10(aves),'o',label='old')
        #ax.errorbar(bin_centers,np.log10(aves1),np.log10(errs),linestyle='None',marker='o',label='new')

        #ax.plot(bin_centers,(aves),'o',label='old')

        #ax.plot(bin_centers,(aves-aves1),'o',label='diff')
        #ax.errorbar(bin_centers,(aves1),(errs),linestyle='None',marker='o',label='new')
        ax.plot(bin_centers,np.log10(aves1),linestyle='None',marker='o',label='new')
        #ax.plot(bin_centers,np.log10(aves1),'o',label='new')
        ax.set_xlabel(r'$b$')
        ax.set_ylabel(r'$\left\langle N_{RC}\right\rangle$')
        ax.set_ylim(0.6,1.6)
        #ax.set_yscale('log')
        #plt.legend()
        plt.savefig('int_rep_2/ave_n_rct.pdf')
        plt.show()
        

    def read_map(self,filename,path='maps/'):
        self.map_data = np.loadtxt(path+filename,delimiter=',')
        #return map_data

    def plot_point_map(self,map_data,axis=6):
        fig, ax = plt.subplots()
    
        #build color bar

        good_data = map_data[:,axis][~np.isnan(map_data[:,axis])]
        cmap_max = np.max(good_data)
        cmap_min = np.min(good_data)
        cmap_max = np.percentile(good_data,95)
        cmap_min = np.percentile(good_data,5)

    
        norm = matplotlib.colors.Normalize(vmin=cmap_min,vmax=cmap_max)
        cmap = matplotlib.cm.get_cmap('cividis_r')

        sc = ax.scatter(map_data[:,0],map_data[:,1],marker='.',c=map_data[:,axis],cmap=cmap,vmin=cmap_min,vmax=cmap_max)#facecolor=cmap(norm(map_data[:,axis])))
        fig.colorbar(sc,ax=ax)#,vmax=cmap_max,vmin=cmap_min)

        ax.set_aspect('equal')
        plt.show()



    def plot_grid_map(self,map_data,lb=True,func=do_nothing,axis=7,cb_label=r'A(K)'):
        fig, ax = plt.subplots()
    
        #build color bar
        good_data = func(map_data[:,axis][~np.isnan(map_data[:,axis])])
        cmap_max = np.max(good_data)
        cmap_min = np.min(good_data)
        cmap_max = np.percentile(good_data,99) 
        cmap_min = np.percentile(good_data,1)

        #norm = matplotlib.colors.Normalize(vmin=cmap_min,vmax=cmap_max)
        #cmap = matplotlib.cm.get_cmap('viridis_r')
        #sc = ax.scatter(np.linspace(cmap_min,cmap_min,100),np.linspace(cmap_min,cmap_min,100)+1e5,c=map_data[:,axis],cmap=cmap,vmin=cmap_min,vmax=cmap_max)#facecolor=cmap(norm(map_data[:,axis])))

        #fig.colorbar(sc,ax=ax)

        norm = matplotlib.colors.Normalize(vmin=cmap_min,vmax=cmap_max)
        cmap = matplotlib.cm.get_cmap('cividis_r')
        #sc = ax.scatter(map_data[:,0],map_data[:,1],facecolor=cmap(norm(map_data[:,axis])))
        sc = ax.scatter(map_data[:,2]+1e3,map_data[:,3],marker='.',c=func(map_data[:,axis]),cmap=cmap,vmin=cmap_min,vmax=cmap_max)#facecolor=cmap(norm(map_data[:,axis])))
        cbar = fig.colorbar(sc,ax=ax)#,vmax=cmap_max,vmin=cmap_min)
        cbar.set_label(cb_label)
    
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
                ax.add_patch(rect)


            ax.set_xlim(2.8,-2.8)
            ax.set_ylim(-2.8,2.8)
            ax.set_xlabel('$l$')
            ax.set_ylabel('$b$')
            
            ax.set_aspect('equal')
            plt.show()

    #def gen_plot_grid_map(self,values,l=self.map_data[:,2],b=self.map_data[:,3],edge_length=self.map_data[:,4],cb_label=r'$A(K)$'):
    def gen_plot_grid_map(self,values,cb_lims=None,cb_label=r'$A(K)$',plotname=None,plotoutlines=False):

        fig,ax = plt.subplots()
        #ipdb.set_trace()
        new_data = np.zeros((len(self.map_data[:,2]),4))
        new_data[:,0] = self.map_data[:,2]
        new_data[:,1] = self.map_data[:,3]
        new_data[:,2] = self.map_data[:,4]
        new_data[:,3] = values
        #new_data = np.zeros(4,len(l))
        #new_data[:,0] = l
        #new_data[:,1] = b
        #new_data[:,2] = edge_length
        #new_data[:,3] = values

        
        good_data = new_data[:,3][~np.isnan(new_data[:,3])]
        if type(cb_lims) == type(None):
            cmap_max = np.max(good_data)
            cmap_min = np.min(good_data)
            cmap_max = np.percentile(good_data,95) 
            cmap_min = np.percentile(good_data,5)

        else:
            cmap_max = cb_lims[1] 
            cmap_min = cb_lims[0]

 
        norm = matplotlib.colors.Normalize(vmin=cmap_min,vmax=cmap_max)
        cmap = matplotlib.cm.get_cmap('cividis_r')
        #sc = ax.scatter(map_data[:,0],map_data[:,1],facecolor=cmap(norm(map_data[:,axis])))
        sc = ax.scatter(new_data[:,0]+1e3,new_data[:,0],marker='.',c=new_data[:,3],cmap=cmap,vmin=cmap_min,vmax=cmap_max)#facecolor=cmap(norm(map_data[:,axis])))
        cbar = fig.colorbar(sc,ax=ax)#,vmax=cmap_max,vmin=cmap_min)
        cbar.set_label(cb_label)


        for pixel in new_data:
            
            #if np.abs(pixel[3]) <1 and pixel[axis]<13.5:
            #    print(pixel[0],',',pixel[1],pixel[2],pixel[3])
            

            if pixel[0]>160:
                rect = pat.Rectangle((pixel[0]-360,pixel[1]),pixel[2],pixel[2],angle=self.angle,facecolor=cmap(norm(pixel[3])))
            else:
                rect = pat.Rectangle((pixel[0],pixel[1]),pixel[2],pixel[2],angle=self.angle,facecolor=cmap(norm((pixel[3]))))
            ax.add_patch(rect)



        if plotoutlines==True:
            fields = []
            #Ugh - numpy doesn't understand why a blank line might mean something
            #Solution by Martin Evans https://stackoverflow.com/questions/36569827/read-txt-data-separated-by-empty-lines-as-several-numpy-arrays
            with open('layout_7f_3.outline.lbad') as fp:
                for k, g in groupby(fp, lambda x: x.startswith(' ')):
                    if not k:
                        fields.append(np.array([[float(x) for x in d.split()] for d in g if len(d.strip())]))
                        
            for f in fields:
                ax.plot(f[:,1],f[:,2],'k-',lw=3)


        ax.set_xlim(2.8,-2.2)
        ax.set_ylim(-2.8,2.1)
        ax.set_xlabel('$l$')
        ax.set_ylabel('$b$')
        
        ax.set_aspect('equal')
        if type(plotname)!=type(None):
            plt.savefig('./int_rep_2/'+plotname)
        plt.show()




    def plot_grid_map(self,lb=True,func=do_nothing,axis=7,cb_label=r'A(K)'):
        fig, ax = plt.subplots()
    
        #build color bar
        good_data = func(self.map_data[:,axis][~np.isnan(self.map_data[:,axis])])
        cmap_max = np.max(good_data)
        cmap_min = np.min(good_data)
        cmap_max = np.percentile(good_data,99) 
        cmap_min = np.percentile(good_data,1)

        #norm = matplotlib.colors.Normalize(vmin=cmap_min,vmax=cmap_max)
        #cmap = matplotlib.cm.get_cmap('viridis_r')
        #sc = ax.scatter(np.linspace(cmap_min,cmap_min,100),np.linspace(cmap_min,cmap_min,100)+1e5,c=map_data[:,axis],cmap=cmap,vmin=cmap_min,vmax=cmap_max)#facecolor=cmap(norm(map_data[:,axis])))

        #fig.colorbar(sc,ax=ax)

        norm = matplotlib.colors.Normalize(vmin=cmap_min,vmax=cmap_max)
        cmap = matplotlib.cm.get_cmap('cividis_r')
        #sc = ax.scatter(map_data[:,0],map_data[:,1],facecolor=cmap(norm(map_data[:,axis])))
        sc = ax.scatter(map_data[:,2]+1e3,map_data[:,3],marker='.',c=func(map_data[:,axis]),cmap=cmap,vmin=cmap_min,vmax=cmap_max)#facecolor=cmap(norm(map_data[:,axis])))
        cbar = fig.colorbar(sc,ax=ax)#,vmax=cmap_max,vmin=cmap_min)
        cbar.set_label(r'$A(K)$')
    
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
                ax.add_patch(rect)


            ax.set_xlim(2.8,-2.8)
            ax.set_ylim(-2.8,2.8)
            ax.set_xlabel('$l$')
            ax.set_ylabel('$b$')
            
            ax.set_aspect('equal')
            plt.show()







if __name__=='__main__':
    #test_map = read_map('maps/test_map1000')
    #tmap = mapper('old_PSF_2017_2.map')
#    tmap = mapper('test2_2017_2.map')
    tmap = mapper('inve_guess_2017_2.map')
    ipdb.set_trace()

    tmap.map_A_K(cb_lims=[0,1.75])
    tmap.map_sigma_H_K()
    #tmap.plot_ave_stars()
    ipdb.set_trace()
    #tmap.plot_ave_stars()
    tmap.plot_N_RC()
    tmap.gen_plot_grid_map(tmap.map_data[:,17])
    tmap.gen_plot_grid_map(tmap.map_data[:,16]) 
    tmap.plot_R()


    tmap.gen_plot_grid_map(tmap.map_data[:,10])

    while True:
        ipdb.set_trace()
        plot_grid_map(test_map,func=do_nothing,axis=7)#[:1000])
        plot_grid_map(test_map,func=A_K,axis=7)#[:1000])
