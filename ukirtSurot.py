#-------------------------------------------------------
# File name: ukirtSurot.py
# 
# Attempts to organize fields by ccd location to create
# larger scale cmds for completeness correction.
#
# Programmer: Aiden Zelakiewicz (zelakiewicz.1@osu.edu)
#
# For questions, reach out to:
#   Aiden Zelakiewicz   (zelakiewicz.1@osu.edu)
#   Samson Johnson      (johnson.7080@osu.edu)
#
# Revision History:
#   21-Mar-2022 :   File created.
#   28-Mar-2022 :   Fixed double correcting the angle.
#   04-Apr-2022 :   Wrapped code into class
#-------------------------------------------------------

#Importing all necessary packages
import os
import pickle
import ipdb
import matplotlib.pyplot as plt
import matplotlib.patches as pat
from astropy import coordinates as coord
from astropy import units as u
import numpy as np
import parameters as pram

class cmd:
    def __init__(self, ra, dec, l=None, b=None, year = pram.year, edge_length=0.25) -> None:

        self.ra = ra
        self.dec = dec

        #Checks if l and b are given. If not, calculate them
        if type(l)==type(None) or type(b)==type(None):
            c = coord.SkyCoord(ra=ra*u.degree,dec=dec*u.degree ,frame='icrs')
            l = c.galactic.l.degree
            b = c.galactic.b.degree
        
        self.l = l
        self.b = b

        self.year = year
        self.edge_length = edge_length

        # Obtaining the UKIRT field positions
        self.fieldData = self.readUKIRTfields()

        self.raw_field_inds = self.findFields(year=self.year, plot_fields=False)

        pass

    def readUKIRTfields(self, filename='ukirtFieldLocations.sav',
                        dir='fieldLocations/'):

        if not os.access(dir+filename, os.R_OK):
            exit('ERROR: no file with CCD locations')

        infile = open(dir+filename,'rb')
        fieldData = pickle.load(infile,encoding='latin1')

        return fieldData
    #Probably not needed anymore...
    def calc_angle(self, ramin,decmin,ramax,decmax, edgelength = 54.4/60.):

        # Calculating the corners of our box
        ra1 = (ramin+ramax)/2 - edgelength/2
        ra2 = (ramin+ramax)/2 - edgelength/2
        dec1 = (decmin+decmax)/2 - edgelength/2
        dec2 = (decmin+decmax)/2 + edgelength/2

        # Creating coordinate objects
        c1 = coord.SkyCoord(ra=ra1*u.degree,dec=dec1*u.degree,frame='icrs')
        c2 = coord.SkyCoord(ra=ra2*u.degree,dec=dec2*u.degree,frame='icrs')

        # Calculates angle
        angle = np.degrees(np.arctan2(c2.galactic.b.degree-c1.galactic.b.degree,c2.galactic.l.degree-c1.galactic.l.degree))
        return angle
    

    def getCCDlocation(self,field,ccd,year=pram.year):
    
        # Each ccd has fov of 13.6'x13.6', so each field has 54.4'x54.4'
        ipdb.set_trace()
        for pos in self.fieldData:
            if pos['year']==str(year) and \
                pos['field']==field and \
                pos['ccd']==ccd:
                RAmin = pos['RAmin']
                RAmax = pos['RAmax']
                DECmin = pos['DECmin']
                DECmax = pos['DECmax']
                print(RAmin," ", RAmax," ", DECmin," ", DECmax)
                return RAmin,RAmax,DECmin,DECmax
        exit('ERROR: year/field/ccd not found in CCD location list')

    def declareFields(self, dim, new_ra=None, new_dec=None):
        """
        Function that gives the index of fields to use based on field size.
        This could be replace later on with a more robust and sensible version
        which returns the fields themselves, but for now it is currently just
        being made to adapt to the old CMD code.

        In its current form, it takes the dimension of the field size 
        """

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

        for i in range(len(self.fieldData)):

            # Makes sure we only obtain fields of same observing year
            if int(self.fieldData[i]["year"]) == int(self.year):

                # This loop will gather the location of fields that the pixel 
                # lies within.
                return
            return

        return

    def getStars(self,skip_read=False,limit_dict=None):


        # Creates an empty limit dictionary if none given
        if type(limit_dict)==type(None):
                limit_dict = {}

        if not skip_read:
        #for each field found in raw field search
            dict_list = []
            for i in raw_field_inds:

                dict_list.append(np.load('../data/ukirt/2017/psfpickles/altBandCrossCheck_'+pram.phot+self.fieldData[i]['year'] \
                    +'_'+self.fieldData[i]['field']+'_'+self.fieldData[i]['ccd']+'.npy',allow_pickle=True).item())


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

            all_dict = dict_list[0]
            for i in np.arange(len(dict_list)-1)+1:
                for key in all_dict.keys():
                    all_dict[key]=np.append(all_dict[key],dict_list[i][key])

        #ipdb.set_trace()
        good_ra_high = all_dict['RA']<=ra_max 
        good_ra_low = all_dict['RA']>=ra_min
        good_ra_inds = good_ra_low & good_ra_high 
        good_dec_high = all_dict['Dec']<=dec_max 
        good_dec_low = all_dict['Dec']>=dec_min
        good_dec_inds = good_dec_low & good_dec_high 
        #good_dec_inds = all_dict['DEC']<=dec_max and all_dict['DEC']>=dec_min
        good_coord_inds = good_ra_inds & good_dec_inds

        #find which catalog entries have all conditions met
        #we will loop through each key and apply that filter, combining with a logical 'and'
        filt_dict = {}
        good_inds = good_coord_inds
        for key in all_dict.keys():
            if key in limit_dict.keys():
                #ipdb.set_trace()
                low_inds = all_dict[key]>=limit_dict[key][0]
                high_inds = all_dict[key]<=limit_dict[key][1]
                try:
                    nan_inds = not (np.isnan(all_dict[key]))
                except:
                    nan_inds = np.ones(len(all_dict[key]),dtype=bool)
                    pass
                good_inds = good_inds & low_inds & high_inds & nan_inds
            else:
                try:
                    good_inds = good_inds & (not np.isnan(all_dict[key]))
                except:
                    nan_inds = np.ones(len(all_dict[key]),dtype=bool)
                    good_inds = good_inds & nan_inds
                    pass
        for key in all_dict.keys():
            filt_dict[key] = all_dict[key][good_inds]
        allStarDict = all_dict
        filterStarDict = filt_dict


        return allStarDict, filterStarDict


if __name__ == '__main__':
    rmin, rmax, dmin, dmax = cmd.getCCDlocation('s5_2', '3')
    print((rmax-rmin))
    print(dmax-dmin)

    fig, ax = plt.subplots()
    ax.plot(rmin,dmin,'.')
    angle = cmd.calc_angle(rmin,dmin,rmax,dmax)
    #rect = pat.Rectangle((rmin,dmin),rmax-rmin,dmax-dmin, angle=-29.97272723566587, facecolor='none',linewidth=1,edgecolor='C0')
    rect = pat.Rectangle((rmin,dmin),rmax-rmin,dmax-dmin, facecolor='none',linewidth=1,edgecolor='C0')
    ax.add_patch(rect)
    ax.set_aspect("equal")
    plt.show()
