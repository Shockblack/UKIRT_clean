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
    def __init__(self, ra, dec, l=None, b=None, year = pram.year, edge_length=0.25):

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

                """
                This loop will gather the location of fields that the pixel 
                lies within. It checks that its location is within the min
                and max of each subfield, then append the index of said field
                into a list of indices for field locations.
                """



                return
            return

        return

    def getStars(self,skip_read=False,limit_dict=None):
        """Function which creates dictionaries of stars which fall under
        the given field ranges. Returns one dictionary which is all the 
        stars in the given field then one for all the stars that pass the
        filter checks.

        Parameters
        ----------
        skip_read : bool, optional
            If True, skips reading in data and checking for correct primary bands.
            Defaults to False
        limit_dict : dict, optional
            A dictionary containing filters on any key found within the stellar
            data. Keys of limit_dict must match keys found in the data files
            and provide lower and upper bound ranges. Primary filters include on 
            number of stars or the magnitude of stars.
            Ex: limit_dict = {'altmag':[12.16],'N':[10,1e6]}
            Defaults to None

        Returns
        -------
        allStarDict : dict
            Dictionary of all stars within the field range, no filters applied.
        filterStarDict : dict
            Dictionary of all stars which pass the filter checks and within the field range.
        """

        # If no limit filter dictionary is given, creates an empty one
        if type(limit_dict)==type(None):
                limit_dict = {}

        if not skip_read:
        
            # Creates an empty list which will hold all the dictionaries of stellar data
            # of the fields found in self.raw_field_inds. This is not the list
            # of stars to be used in the fit, but every star in the given fields.
            dict_list = []

            # Loops over each field found to contain stars in the ra and dec range
            # given by the findfields() function. Stored in self.raw_field_inds
            for i in self.raw_field_inds:

                # Appending the dictionary list with all the stellar data
                dict_list.append(np.load('../data/ukirt/2017/psfpickles/altBandCrossCheck_'+pram.phot+self.fieldData[i]['year'] \
                    +'_'+self.fieldData[i]['field']+'_'+self.fieldData[i]['ccd']+'.npy',allow_pickle=True).item())

                # Below is a check to make sure out primary band is always the H-band and that the secondary is K.
                # That is, we want to check if 'mag' is H and if 'altmag' is K.

            for i in range(len(dict_list)):
                # Creating a dictionary of the given iteration used to compare against
                tdict = dict_list[i]

                # Checking if the 'delta' key is not present and adding it in if that's the case.
                if 'delta' not in tdict.keys():
                    tdict['delta'] = tdict['mag'] - tdict['altmag']
                

                # If the name of the first entry has _K_, then K is the primary band.
                # We want to switch these if that is the case so that H is primary.
                if '_K_' in tdict[list(tdict.keys())[0]][0]:

                    # Creates a temporary dictionary for restructuring
                    tempdict = {}

                    for key in tdict.keys():

                        # If K is primary band, then delta (K-H) is the negative of the delta (H-K) which we want.
                        if key == 'delta':
                            tempdict[key] = -tdict[key]

                        # These are fine as is, so just copy over.
                        elif key == 'error' or key == 'offset':
                            tempdict[key]=tdict[key]

                        # If the prefix "alt" is not in the key then it corresponds to the K band.
                        # We want to flip this so that K IS the alt band by adding "alt" before each key.
                        elif 'alt' not in key:
                            tempdict['alt'+key]=tdict[key]

                        # This means that everything left is 'alt'+key, which here corresponds to the H-band values.
                        # We want this to be just key, without the 'alt.' To do this we grab all the remaining keys
                        # and splice the string name so that the first 3 letters 'alt' get removed.
                        else:
                            tempdict[key[3:]]=tdict[key]
                
                        # Set the temporary dictionary dict as the entry
                        dict_list[i]=tempdict

            # Here the dictionary lists are restructured. Instead of being a list of dictionaries,
            # this will turn it into a dictionary of lists. This line below creates the "framework"
            # for the dictionary by pulling the first entry of the dict_list[0] so that the keys
            # are already in place.
            all_dict = dict_list[0]

            # Creates a range to loop over all the stars excluding the first one. The -1 in the arange
            # after the len() makes it 1 item shorter. The +1 after the arange adds 1 to every
            # value in the arange, essentially making it start 1 index later.
            for i in np.arange(len(dict_list)-1)+1:
                # Loops over all keys
                for key in all_dict.keys():
                    # Appends the dictionary key with the value of the next datapoint of the i-th star
                    all_dict[key]=np.append(all_dict[key],dict_list[i][key])

        # Gathers the indices for all stars within the RA and DEC range.
        good_ra_high = all_dict['RA']<=self.ra_max 
        good_ra_low = all_dict['RA']>=self.ra_min
        good_ra_inds = good_ra_low & good_ra_high

        good_dec_high = all_dict['Dec']<=self.dec_max 
        good_dec_low = all_dict['Dec']>=self.dec_min
        good_dec_inds = good_dec_low & good_dec_high 
        
        good_coord_inds = good_ra_inds & good_dec_inds

        # Need to find entries which pass all conditions. Will loop through each key, applying
        # the filter and adjusting the "good indices" for each.

        # Creating empty filtered star dictionary
        filt_dict = {} 
        # Starting off with our good indices being those within coord range
        good_inds = good_coord_inds
        # Looping over all keys
        for key in all_dict.keys():

            # Checks if the specified key is within our limiting filter dictionary.
            if key in limit_dict.keys():

                # Gathers indices where the stars in the dictionary pass the filter.
                low_inds = all_dict[key]>=limit_dict[key][0]
                high_inds = all_dict[key]<=limit_dict[key][1]

                # Gathers all elements that don't have "nan" entries
                try:
                    # Returns False if 'nan' and True if not, basically opposite of Numpy.isnan
                    nan_inds = not (np.isnan(all_dict[key]))
                except:
                    # If there are no 'nan' values, return a list of all True
                    nan_inds = np.ones(len(all_dict[key]),dtype=bool)
                    pass

                # Uses logical 'and' to gather all indices that have passed checks
                good_inds = good_inds & low_inds & high_inds & nan_inds

            # If no filter key, just checks for 'nan' same way as above
            else:
                try:
                    good_inds = good_inds & (not np.isnan(all_dict[key]))
                except:
                    nan_inds = np.ones(len(all_dict[key]),dtype=bool)
                    good_inds = good_inds & nan_inds
                    pass
        
        # Appends our filtered dictionary with values in 'good_inds'
        for key in all_dict.keys():
            filt_dict[key] = all_dict[key][good_inds]
        allStarDict = all_dict
        filterStarDict = filt_dict

        # Returns the dictionaries
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
