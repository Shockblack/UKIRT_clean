from operator import le
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as pat
#from matplotlib import cm
from collections import defaultdict
import ipdb
from readData import readUKIRTfields
from readData import readAltbandCrossCheck
import astropy
from astropy import coordinates as coord
from astropy import units as u
import parameters as pram
import RedClumpFinder as rcf

class cmd:
    def __init__(self,filename,ra,dec,l=None,b=None,edge_length=0.25,year=pram.year,just_fields=False):
        #set fundamental box properties, center and length of box
        self.ra = ra
        self.dec = dec
        if type(l)==type(None) or type(b)==type(None):
            c = coord.SkyCoord(ra=ra*u.degree,dec=dec*u.degree ,frame='icrs')
            l = c.galactic.l.degree
            b = c.galactic.b.degree
        self.l = l
        self.b = b

        self.edge_length = edge_length
        self.year=year

        #set the upper and lower limits for our box in ra and dec
        self.ra_max = self.ra+self.edge_length/2.
        self.ra_min = self.ra-self.edge_length/2.
        self.dec_max = self.dec+self.edge_length/2.
        self.dec_min = self.dec-self.edge_length/2.
        
        #read in the field data
        self.fieldData = readUKIRTfields()        
        

        self.raw_field_inds = self.findFields(year=2017,plot_fields=False)

        if just_fields: 
            return


        self.limit_dict = {'N':[10,1e6],'altN':[3,1e6],'offset':[-0.1,0.1], 'altMAD':[-0.1,0.1], 'MAD':[-0.1,0.1]}# {'mag':[12,15]}
        self.getStars(limit_dict=self.limit_dict)
        self.cm_dict = {'altmag':[12,16],'delta':[-1,5]}
        #self.cm_dict = {'altmag':[12,16],'delta':[-1,5]}
        self.color_mag_cut(self.cm_dict,percentile_cut=True)
        self.histMag(plot=False)
        #ra_col = 1
        #dec_col = 2
        #self.raw_data = np.genfromtxt(filename,skip_header=1,delimiter='|')
        #self.data_ra_low = self.raw_data[self.raw_data[:,ra_col]]

    
    def fetchStars(self,ra,dec,edge_length=pram.arcmin/60.,year=pram.year):
        return 0
    
    def updateCMD(self,new_ra,new_dec,l=None,b=None):
        #set fundamental box properties, center and length of box
        self.ra = new_ra
        self.dec = new_dec
        if type(l)==type(None) or type(b)==type(None):
            c = coord.SkyCoord(ra=new_ra*u.degree,dec=new_dec*u.degree ,frame='icrs')
            l = c.galactic.l.degree
            b = c.galactic.b.degree
        self.l = l
        self.b = b
        #self.edge_length = edge_length
        #self.year=year
        
        #set the upper and lower limits for our box in ra and dec
        self.ra_max = self.ra+self.edge_length/2.
        self.ra_min = self.ra-self.edge_length/2.
        self.dec_max = self.dec+self.edge_length/2.
        self.dec_min = self.dec-self.edge_length/2.

        #read in the field data
        #self.fieldData = readUKIRTfields()
        
        self.raw_field_inds = self.findFields(year=2017)
        #self.limit_dict = {'N':[10,1e6],'altN':[3,1000]}# {'mag':[12,15]}
        self.getStars(skip_read=True,limit_dict=self.limit_dict)
        #self.cm_dict = {'altmag':[12,16],'delta':[-1,5]}
        self.color_mag_cut(self.cm_dict,percentile_cut=True)
        self.histMag(plot=False)
        
        print('CMD updated!')

    def check_new_coords(self,new_ra,new_dec):
        new_raw_fields = self.findFields(new_ra=new_ra,new_dec=new_dec,year=self.year)
        if np.array_equal(np.array(new_raw_fields),np.array(self.raw_field_inds)):
            return True
        else:
            return False


    def findFields(self,year=pram.year,new_ra=None,new_dec=None,plot_fields=True):#,ra=self.ra,dec=self.dec,edge_length=0.25):
        #logic for checking if we need to reload
        if new_ra == None and new_dec == None:
            ra_max = self.ra_max
            ra_min = self.ra_min
            dec_max = self.dec_max
            dec_min = self.dec_min
        else:
            ra_max = new_ra + self.edge_length/2.
            ra_min = new_ra - self.edge_length/2.
            dec_max = new_dec + self.edge_length/2.
            dec_min = new_dec - self.edge_length/2.
        #determine which fields our cmd will use
        raw_field_inds = []
        for i in range(len(self.fieldData)):
            #check to see if included in year we want. 
            if int(self.fieldData[i]['year'])==int(year):
                #ipdb.set_trace()
                #check if the field lies in the RA range
                #if (self.fieldData[i]['RAmax']>=self.ra_min and self.fieldData[i]['RAmax']<=self.ra_max) or \
                #        (self.fieldData[i]['RAmin']>=self.ra_min and self.fieldData[i]['RAmin']<=self.ra_max):

                ###NEED TO MAKE THIS BASED ON MINIFIELD LIMIT
                if (ra_min<=self.fieldData[i]['RAmax'] and ra_min>=self.fieldData[i]['RAmin']) or \
                        (ra_max<=self.fieldData[i]['RAmax'] and ra_max>=self.fieldData[i]['RAmin']) or\
                        (ra_max>=self.fieldData[i]['RAmax'] and ra_min<=self.fieldData[i]['RAmin']):

                    #check if the field lies in the DEC range
                    if (dec_min<=self.fieldData[i]['DECmax'] and dec_min>=self.fieldData[i]['DECmin']) or \
                            (dec_max<=self.fieldData[i]['DECmax'] and dec_max>=self.fieldData[i]['DECmin']) or \
                            (dec_max>=self.fieldData[i]['DECmax'] and dec_min<=self.fieldData[i]['DECmin']):
                    #if (self.fieldData[i]['DECmax']>=self.dec_min and self.fieldData[i]['DECmax']<=self.dec_max) or \
                    #        (self.fieldData[i]['DECmin']>=self.dec_min and self.fieldData[i]['DECmin']<=self.dec_max):
                        #append the fieldData ind i to our list
                        raw_field_inds.append(i)
                        #break
        #self.raw_field_inds =  raw_field_inds
        return raw_field_inds
        #return 'Found raw field inds'



    def getStars(self,skip_read=False,limit_dict=None):
        """
        filter_dict must have key match those read in form files, and must provide lower and higher bound ranges. 
        e.g. for magnitude it will be filter_dict['mag']=[low_mag, high_mag]
        """
        #catch if there is no filter dictionary
        
        if type(limit_dict)==type(None):
            limit_dict = {}
        #ipdb.set_trace()
        if not skip_read:
        #for each field found in raw field search
            dict_list = []
            for i in self.raw_field_inds:

                #SAMSON: if you are trying to run, you need to remove everything before psfpickles -Aiden
                #I also added an allow_pickle argument to be True
                #added photom type from parameter file... If broken just remove pram.phot with 'PSF'
                dict_list.append(np.load('../data/ukirt/2017/psfpickles/altBandCrossCheck_'+pram.phot+self.fieldData[i]['year']+'_'+self.fieldData[i]['field']+'_'+self.fieldData[i]['ccd']+'.npy',allow_pickle=True).item())
                

        
                #going to put a check to make sure we are consistent with which 'mag' and 'altmag' refer to H and K. 
                #H is always primary, K is always alt
                #this is a bit cludgy, but couldn't think of a more graceful way 
            for i in range(len(dict_list)):
                tdict = dict_list[i]
                #ipdb.set_trace()
                #need to add catch if 'delta' key is not present
                if 'delta' not in tdict.keys():
                    tdict['delta'] = tdict['mag'] - tdict['altmag']

                #if the name of the first entry has _K_, then K is the primary band.
                #we want to switch these
                if '_K_' in tdict[list(tdict.keys())[0]][0]:
                    tempdict = {}
                    for key in tdict.keys():
                    #the delta (K-H) is the negative of what we want(H-K), so fix that
                        if key == 'delta':
                            tempdict[key] = -tdict[key]
                    #these are fine as is, so just copy
                        elif key == 'error' or key == 'offset':
                            tempdict[key]=tdict[key]
                    #if 'alt' is not in the key it corresponds to K, so we want that as 'alt'+key
                        elif 'alt' not in key:
                            tempdict['alt'+key]=tdict[key]
                    #otherwise, 'alt'+key is the H-band value and we want that as just key
                    #here, 'altname'[3:] is just 'name', the slice removes the preceeding 'alt'
                        else:
                            tempdict[key[3:]]=tdict[key]
                #ipdb.set_trace()
                #put this tempdict at the entry
                        dict_list[i]=tempdict

                

            self.all_dict = dict_list[0]
            for i in np.arange(len(dict_list)-1)+1:
                for key in self.all_dict.keys():
                    self.all_dict[key]=np.append(self.all_dict[key],dict_list[i][key])

        

        #ipdb.set_trace()
        good_ra_high = self.all_dict['RA']<=self.ra_max 
        good_ra_low = self.all_dict['RA']>=self.ra_min
        good_ra_inds = good_ra_low & good_ra_high 

        good_dec_high = self.all_dict['Dec']<=self.dec_max 
        good_dec_low = self.all_dict['Dec']>=self.dec_min
        good_dec_inds = good_dec_low & good_dec_high 

        #good_dec_inds = all_dict['DEC']<=self.dec_max and all_dict['DEC']>=self.dec_min
        good_coord_inds = good_ra_inds & good_dec_inds
        
        #find which catalog entries have all conditions met
        #we will loop through each key and apply that filter, combining with a logical 'and'
        filt_dict = {}
        good_inds = good_coord_inds
        for key in self.all_dict.keys():
            if key in limit_dict.keys():
                #ipdb.set_trace()
                low_inds = self.all_dict[key]>=limit_dict[key][0]
                high_inds = self.all_dict[key]<=limit_dict[key][1]
                try:
                    nan_inds = not (np.isnan(self.all_dict[key]))
                except:
                    nan_inds = np.ones(len(self.all_dict[key]),dtype=bool)
                    pass
                good_inds = good_inds & low_inds & high_inds & nan_inds
            else:
                try:
                    good_inds = good_inds & (not np.isnan(self.all_dict[key]))
                except:
                    nan_inds = np.ones(len(self.all_dict[key]),dtype=bool)
                    good_inds = good_inds & nan_inds
                    pass


        for key in self.all_dict.keys():
            filt_dict[key] = self.all_dict[key][good_inds]

        self.allStarDict = self.all_dict
        self.filterStarDict = filt_dict

        #XXX still need a way to cross check stars, make sure overlapping fields don't duplicate stars
        #names are not unique. Will most likelty need to look for RA,DEC distance less than some threshold. 
        return

    def color_mag_cut(self,limit_dict=None,percentile_cut=False):
        if type(limit_dict)==type(None):
            print('No limits on color or magnitude')

        #make a boolean index array to initialize our color_mag cuts
        good_inds = np.ones(len(self.filterStarDict['mag']),dtype=bool) 

        #apply magnitude limit
        if 'mag' in limit_dict.keys():
            low_inds = self.filterStarDict['mag']>=limit_dict['mag'][0]
            high_inds = self.filterStarDict['mag']<=limit_dict['mag'][1]
            good_inds = low_inds & high_inds

        #apply alt-magnitude limit
        if 'altmag' in limit_dict.keys():
            low_inds = self.filterStarDict['altmag']>=limit_dict['altmag'][0]
            high_inds = self.filterStarDict['altmag']<=limit_dict['altmag'][1]
            good_inds = good_inds & low_inds & high_inds
        #ipdb.set_trace()
        if 'delta' in limit_dict.keys():
            if percentile_cut:
                limit_dict['delta'][0]=np.percentile(self.filterStarDict['delta'],15) #only including stars above 15th percentile in (H-K) color
            low_inds = (self.filterStarDict['delta']) >= limit_dict['delta'][0]
            high_inds = (self.filterStarDict['delta']) <= limit_dict['delta'][1]

            good_inds = good_inds & low_inds & high_inds
        
        self.fitStarDict = {}
        for key in self.filterStarDict.keys():
            self.fitStarDict[key] = self.filterStarDict[key][good_inds]
        return 'Stars cut by magnitude and color'


    def color_fit_cut(self,M_RC, sigma_RC,Nsigma=2):
        good_inds = np.ones(len(self.filterStarDict['mag']),dtype=bool)
        
        #cut by percentile in color
        low_inds = (self.filterStarDict['delta']) >=np.percentile(self.filterStarDict['delta'],20)
        good_inds = good_inds & low_inds

        #only grab stars within some mulitple of sigma of the red clump
        low_lim = M_RC-Nsigma*sigma_RC
        high_lim = M_RC+Nsigma*sigma_RC
        low_inds = (self.filterStarDict['altmag']) >= low_lim
        high_inds = (self.filterStarDict['altmag']) <= high_lim
        good_inds = good_inds & low_inds & high_inds

        self.colorFitStarDict = {}
        for key in self.filterStarDict.keys():
            self.colorFitStarDict[key] = self.filterStarDict[key][good_inds]
        return 'Stars cut by magnitude and color'

    def plotFields(self):
        fig, ax = plt.subplots()
        ax.plot(self.ra_min,self.dec_min,'.')
        rect = pat.Rectangle((self.ra_min,self.dec_min),self.edge_length,self.edge_length,facecolor='none',linewidth=1,edgecolor='C1')
        ax.add_patch(rect)
        for i in self.raw_field_inds:
            print(i,self.fieldData[i]['RAmin'])
            ax.plot(self.fieldData[i]['RAmin'],self.fieldData[i]['DECmin'],'.')
            rect = pat.Rectangle((self.fieldData[i]['RAmin'],self.fieldData[i]['DECmin']),self.fieldData[i]['RAmax']-self.fieldData[i]['RAmin'], \
                                     self.fieldData[i]['DECmax']-self.fieldData[i]['DECmin'], facecolor='none',linewidth=1,edgecolor='C0')
            ax.add_patch(rect)
        try:
            #ipdb.set_trace()
            ax.plot(self.filterStarDict['RA'],self.filterStarDict['Dec'],'.',markersize=1)
        except:
            print('No stars to plot!\n')

        ax.set_xlabel('Right Ascension')
        ax.set_ylabel('Declination')
        plt.show()


    def plotCMD(self,cm_dict=None,fitcorr=True):
        #plt.rc('text',usetex=True)
        fig,ax = plt.subplots()
        c = coord.SkyCoord(ra=self.ra*u.degree,dec=self.dec*u.degree,frame='icrs')

        ax.set_title(r'($l,b$)=(%.2f,%.2f), size=%0.1f$^\prime$'%(c.galactic.l.degree,c.galactic.b.degree,60*self.edge_length))

        ax.plot(self.filterStarDict['delta'],self.filterStarDict['altmag'],'.',markersize=2)
        try:

            ax.plot(self.fitStarDict['delta'],self.fitStarDict['altmag'],'.',markersize=2)
        except:
            print('No stars for fit dictionary')
        try:

            ax.plot(self.colorFitStarDict['delta'],self.colorFitStarDict['altmag'],'.',markersize=2)
            if fitcorr:
                ipdb.set_trace()
                mean_delta = np.mean(self.colorFitStarDict['delta'])
                mean_altmag = np.mean(self.colorFitStarDict['altmag'])

                sigma_delta = np.std(self.colorFitStarDict['delta'])
                sigma_altmag = np.std(self.colorFitStarDict['altmag'])
                
                corr = np.sum((1./len(self.colorFitStarDict['delta']))*(self.colorFitStarDict['delta']-mean_delta)*(self.colorFitStarDict['altmag']-mean_altmag))/sigma_delta/sigma_altmag
                print('Correlation: %f\n'%corr)
                
                coeffs = np.polyfit(self.colorFitStarDict['delta'],self.colorFitStarDict['altmag'],1)
                
                fitx = np.linspace(self.colorFitStarDict['delta'].min(),self.colorFitStarDict['delta'].max(),100)
                ax.plot(fitx,np.polyval(coeffs,fitx))

        except:
            print('No stars for color fit dictionary')
        if type(cm_dict)!=type(None):
            limit_box = []
            #ipdb.set_trace()
            keys = []
            for key in cm_dict.keys():
                keys.append(key)
            limit_box.append([cm_dict[keys[0]][0],cm_dict[keys[0]][1],cm_dict[keys[0]][1],cm_dict[keys[0]][0],cm_dict[keys[0]][0]])
            limit_box.append([cm_dict[keys[1]][0],cm_dict[keys[1]][0],cm_dict[keys[1]][1],cm_dict[keys[1]][1],cm_dict[keys[1]][0]])
            for i in range(len(limit_box[0])):
                print(limit_box[1][i],limit_box[0][i])
            ax.plot(limit_box[1],limit_box[0],'C1',linewidth=2)
            #for key in cm_dict.keys():

        #plot some vertical percentile lines
        #for perc in np.linspace(20,80,7):
        #    ax.axvline(np.percentile(self.filterStarDict['delta'],perc))

        ax.set_xlabel(r'$\it{H-K}$')
        ax.set_ylabel(r'$\it{K}$')
        lower,upper = np.percentile(self.filterStarDict['delta'],[5,95])
        #ax.set_xlim([lower,upper])
        ax.set_xlim([-.3,2])
        ax.set_ylim([18,11])
        #ax.axvline(0.09)
        #plt.gca().invert_yaxis()
        #plt.savefig('../misc_figs/cmd1.pdf')
        plt.show()

    def plotCMDhist(self,cm_dict=None,fitcorr=True, fit=None, plotsave=False,figdir='../misc_figs/CMDSideHist.pdf',pixID=None):
        #Function plots a CMD with a histogram subplot using gridspec
        from matplotlib.gridspec import GridSpec
        gs = GridSpec(7, 3)

        fig = plt.figure()
        ax = fig.add_subplot(gs[0:7, 0:2])

        ax.set_xlabel(r'$\it{H-K}$')
        ax.set_ylabel(r'$\it{K}$')
        #ipdb.set_trace()
        percmax = np.percentile(self.filterStarDict['delta'],99)
        ax.set_xlim([-.3,percmax]) #TEST THIS
        ax.set_ylim([18,11])

        ax.plot(self.filterStarDict['delta'],self.filterStarDict['altmag'],'.',markersize=2,color='silver')
        try:
            #ax.plot(self.fitStarDict['mag']-self.fitStarDict['altmag'],self.fitStarDict['altmag'],'.',markersize=1)
            #Data enclosed in black area, used for luminosity fit
            ax.plot(self.fitStarDict['delta'],self.fitStarDict['altmag'],'.',markersize=2,color='black')
        except:
            print('No stars for fit dictionary')

        numbins = max([20,int(5*np.log(len(self.fitStarDict['altmag'])))])
        #Code for the side figure
        ax_hist_y = fig.add_subplot(gs[1:5, 2])
        ax_hist_y.set_ylim([16,12])

        try:
            #ax.plot(self.fitStarDict['mag']-self.fitStarDict['altmag'],self.fitStarDict['altmag'],'.',markersize=1)
            ax.plot(self.colorFitStarDict['delta'],self.colorFitStarDict['altmag'],'.',markersize=2)
            if fitcorr:
                ipdb.set_trace()
                mean_delta = np.mean(self.colorFitStarDict['delta'])
                mean_altmag = np.mean(self.colorFitStarDict['altmag'])

                sigma_delta = np.std(self.colorFitStarDict['delta'])
                sigma_altmag = np.std(self.colorFitStarDict['altmag'])
                
                corr = np.sum((1./len(self.colorFitStarDict['delta']))*(self.colorFitStarDict['delta']-mean_delta)*(self.colorFitStarDict['altmag']-mean_altmag))/sigma_delta/sigma_altmag
                print('Correlation: %f\n'%corr)
                
                coeffs = np.polyfit(self.colorFitStarDict['delta'],self.colorFitStarDict['altmag'],1)
                
                fitx = np.linspace(self.colorFitStarDict['delta'].min(),self.colorFitStarDict['delta'].max(),100)
                ax.plot(fitx,np.polyval(coeffs,fitx),color='black')

        except:
            print('No stars for color fit dictionary')

        
        if type(cm_dict)!=type(None):
            limit_box = []
            #ipdb.set_trace()
            keys = []
            for key in cm_dict.keys():
                keys.append(key)
            limit_box.append([cm_dict[keys[0]][0],cm_dict[keys[0]][1],cm_dict[keys[0]][1],cm_dict[keys[0]][0],cm_dict[keys[0]][0]])
            limit_box.append([cm_dict[keys[1]][0],cm_dict[keys[1]][0],cm_dict[keys[1]][1],cm_dict[keys[1]][1],cm_dict[keys[1]][0]])
            for i in range(len(limit_box[0])):
                print(limit_box[1][i],limit_box[0][i])
            ax.plot(limit_box[1],limit_box[0],color='black',linewidth=2)

        fitinds = self.fitMagHist != 0
        yvals = self.fitMagHist[fitinds]
        xvals = ((self.magBins[1:]+self.magBins[:-1])/2)[fitinds]

        ax_hist_y.barh(xvals,yvals,color='dimgray',height=0.5*4/numbins)
        if type(fit) != type(None):
            ax_hist_y.plot(fit,xvals,color='black',linewidth=2)
            pass
        ax_hist_y.tick_params(axis="y", labelleft=False)
        
        #Sets limits on the histogram by utilizing rounding
        cap = np.round(np.array(yvals.max(),dtype=float) / 50) *50
        if cap < yvals.max():
            cap += 50
        #ax_hist_y.set_xticks(np.arange(0,cap,step=50))
        inval = int(1.1*yvals.max()/4)
        ax_hist_y.set_xticks(np.arange(0,inval*4,inval))

        if self.l > 350:
            self.l -=360

        if type(pixID) != type(None): #For if I want to know know what pixel this specific plot comes from. Useful for debugging
            ax.set_title('(l,b)=(%.4f,%.4f), size=%0.1f$^\prime$, Pix=%d'%(self.l,self.b,pram.arcmin,pixID))
        else:
            ax.set_title('(l,b)=(%.4f,%.4f), size=%0.1f$^\prime$'%(self.l,self.b,pram.arcmin))
        plt.subplots_adjust(wspace=0.05) #Makes plots close together
        if plotsave == True: #For use when running checks. ALlows a code to continue running by just setting it to save
            plt.savefig(figdir)
        elif plotsave == False:
            plt.show()

    def histCMD(self,cm_dict=None,binnumber=100):
        
        lower,upper = np.percentile(self.filterStarDict['delta'],[5,95])
        delta_bins = np.linspace(lower,upper,binnumber)
        k_bins = np.linspace(11,18,binnumber)
        H, xedges, yedges = np.histogram2d(self.filterStarDict['delta'], self.filterStarDict['altmag'], (delta_bins,k_bins))

        #plt.rc('text',usetex=True)
        fig,ax = plt.subplots()
        p = ax.pcolormesh(xedges,yedges,((H).T))
        ax.set_aspect('auto')
        #cbar = plt.colorbar(p)
        levels = np.array([0.05,0.25,0.5,0.75,0.95])*H.max()
        #ax.contour(H.T,extent=[xedges.min(),xedges.max(),yedges.min(),yedges.max()],levels=levels,linewidths=3,cmap=plt.get_cmap('viridis_r'))
        for perc in np.linspace(20,80,7):
            ax.axvline(np.percentile(self.filterStarDict['delta'],perc),color='C1')
        #fig2 = plt.figure()
        

        #plt.hist2d(self.filterStarDict['delta'], self.filterStarDict['altmag'], bins=100)

        
        plt.xlabel('H-K')
        plt.ylabel('K')
        #cbar = plt.colorbar()
        #cbar.ax.set_ylabel('Counts')

        plt.gca().invert_yaxis()
        #plt.show()

        #ipdb.set_trace()
        return 0
        #ax.plot(self.filterStarDict['mag']-self.filterStarDict['altmag'],self.filterStarDict['altmag'],'.',markersize=1)
        ax.plot(self.filterStarDict['delta'],self.filterStarDict['altmag'],'.',markersize=1)
        try:
            #ax.plot(self.fitStarDict['mag']-self.fitStarDict['altmag'],self.fitStarDict['altmag'],'.',markersize=1)
            ax.plot(self.fitStarDict['delta'],self.fitStarDict['altmag'],'.',markersize=1)
        except:
            print('No stars for fit dictionary')
        if type(cm_dict)!=type(None):
            limit_box = []
            #ipdb.set_trace()
            keys = []
            for key in cm_dict.keys():
                keys.append(key)
            limit_box.append([cm_dict[keys[0]][0],cm_dict[keys[0]][1],cm_dict[keys[0]][1],cm_dict[keys[0]][0],cm_dict[keys[0]][0]])
            limit_box.append([cm_dict[keys[1]][0],cm_dict[keys[1]][0],cm_dict[keys[1]][1],cm_dict[keys[1]][1],cm_dict[keys[1]][0]])
            for i in range(len(limit_box[0])):
                print(limit_box[1][i],limit_box[0][i])
            ax.plot(limit_box[1],limit_box[0],'C1')
            #for key in cm_dict.keys():

        ax.set_xlabel(r'$\it{H-K}$')
        ax.set_ylabel(r'$\it{K}$')
        ax.set_xlim([-.2,1])
        ax.set_ylim([19,10])
        ax.axvline(0.09)
        #plt.gca().invert_yaxis()
        plt.show()

    def histMag(self,binnumber=100,plot=False):
        #ipdb.set_trace()
        bins = np.linspace(11,18,binnumber)
        self.totMagHist, self.magBins = np.histogram(self.filterStarDict['altmag'],bins)
        
        #Might want to comment out next two lines if overlapping occurs
        vals = self.totMagHist
        bins = self.magBins

        #vals, bins = np.histogram(self.filterStarDict['altmag'],bins)
        bin_centers = (bins[1:]+bins[:-1])/2
        

        try:
            #ipdb.set_trace()
            self.fitMagHist,bins2 = np.histogram(self.fitStarDict['altmag'],bins)
            #Might want to comment out next two lines if overlapping occurs

            vals1 = self.fitMagHist
            vals2 = vals/vals.max()-vals1/vals.max()
        except:

            print('No stars for fit to plot')

        if plot:
            fig,ax = plt.subplots()
            ax.plot(bin_centers,vals/vals.max())    
            ax.plot(bin_centers,vals1/vals.max())
            ax.plot(bin_centers,vals2)#/vals2.max())
            ax.set_xlabel(r'$\it{K}$')
            plt.show()

    def histColor(self,binnumber=100,plot=True):
        #ipdb.set_trace()
        bins = np.linspace(-.3,1.5,binnumber)
        vals,bins = np.histogram(self.filterStarDict['delta'],bins)
        bin_centers = (bins[1:]+bins[:-1])/2
        
        fig,ax = plt.subplots()
        ax.axvline(0.09)
        #ax.plot(bin_centers,vals/vals.max())
        ax.plot(bin_centers,np.cumsum(vals/vals.max()))
        try:
            vals1,bins1 = np.histogram(self.fitStarDict['delta'],bins)
            #ax.plot(bin_centers,vals1/vals.max())
            ax.plot(bin_centers,np.cumsum(vals1/vals.max()))
            vals2 = vals/vals.max()-vals1/vals.max()
            #ax.plot(bin_centers,np.cumsum(vals2))#/vals2.max())
            ax.plot(bin_centers,np.cumsum(vals2))#/vals2.max()))
        except:
            print('No stars for fit to plot')
        ax.set_xlabel(r'$\it{(H-K)}$')
        ax.set_yscale('log')
        #plt.show()

if __name__=='__main__':
    
    #ipdb.set_trace()
    #cmd_test = cmd('test.txt',268.5,-28.7,edge_length=240/60.)

    # Creating plots for the paper using these two locations
    # Center field
    cmd_test = cmd('test.txt',266.,-29.,edge_length=pram.arcmin/60.)
    # Edge field
    # cmd_test = cmd('test.txt',269.,-29.5,edge_length=pram.arcmin/60.)

    #cmd_test = cmd('test.txt',269.36695 , -28.98901,edge_length=pram.arcmin/60.)
    #cmd_test = cmd('test.txt',268.5,-29.7,edge_length=.1)

    #cmd_test = cmd('test.txt',267.,-30.5,edge_length=.1)
    #cmd_test = cmd('test.txt',269.,-29.3,edge_length=.05)
    #cmd_test = cmd('test.txt',267.,-28.5,edge_length=.05)
    cmd_test.findFields(year=pram.year)
    #cmd_test.plotFields()
    #ipdb.set_trace()
    limit_dict = {'N':[10,1e6],'altN':[3,1000]}# {'mag':[12,15]}
    cmd_test.getStars(limit_dict)

    # import csv
    # with open('../data/cmd_example.csv', 'w') as file:
    #     # 2. step
    #     writer = csv.writer(file)
    #     # 3. step
    #     for i in range(len(cmd_test.filterStarDict['altmag'])):
    #         writer.writerow([cmd_test.filterStarDict['altmag'][i],cmd_test.filterStarDict['delta'][i]])

    
    cm_dict = {'altmag':[12,16],'delta':[-1,5]}
    cmd_test.color_mag_cut(cm_dict,percentile_cut=True)

    #ipdb.set_trace()
    #cmd_test.plotFields()
    #plt.savefig('plots_sam/fieldplot_sample.2.pdf')
    #plt.show()
    rcfinder=rcf.redclumpfinder(cmd_test)
    M_RCguess = rcfinder.icMethod()
    fit = rcfinder.fitRCMagnitude(rcf.redclumpOnlyExp, plotfit=False, M_RC=M_RCguess)
    cmd_test.plotCMDhist(cm_dict,fit=fit, plotsave=True, figdir='../misc_figs/CMD_Example_High.pdf')
    #cmd_test.histCMD(cm_dict,binnumber=25)
    #cmd_test.histMag(binnumber=100)
    #cmd_test.histColor(binnumber=100)
    #plt.show()

    
    #plt.savefig('plots_sam/CMD_sample.2.pdf')
    #ipdb.set_trace()
    #test = [np.array(u) for u in set([tuple(j,k) for j,k in cmd_test.allStarDict['RA'],cmd_test.allStarDict['Dec']])]
    #data = readAltbandCrossCheck(cmd_test.fieldData[cmd_test.fields[0]]['year'],cmd_test.fieldData[cmd_test.fields[0]]['field'],cmd_test.fieldData[cmd_test.fields[0]]['ccd'])
    
    #I commented these out, don't know why they here. Just printed long list of errors
    #fieldCenters = readUKIRTfields()
    
    #path = '/mnt/a/documents/files/surp/data/ukirt/2017/altBandCrossCheck/' #File path used by Aiden
    #path = '/Users/johnson.7080/work/jplsip/extinctionMaps/altBandCrossCheck/' #Samson file path
    #filename = 'altBandCrossCheck_CASU2018_s6_4_4.txt'

    #raw = np.genfromtxt(path+filename,skip_header=1,delimiter='|')

    

    #ipdb.set_trace()
