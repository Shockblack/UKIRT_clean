#-----------------------------------------------------------------------
# File name: createCMD.py
# 
# A file adapted from original plotCMD.py file which gathers photometry
# data to create CMDs. This also has started organization to derive 
# reddening vectors for our fields. This is in the process in moving away
# from luminosity fitting similar to Nataf et al. 2013 and doing the
# weight method like Surot et al. 2020.
#
# Programmer: Aiden Zelakiewicz (zelakiewicz.1@osu.edu)
#             Samson Johnson    (johnson.7080@osu.edu)
#
# For questions, reach out to:
#   Aiden Zelakiewicz   (zelakiewicz.1@osu.edu)
#   Samson Johnson      (johnson.7080@osu.edu)
#
# Revision History:
#   21-Mar-2022 :   File created.
#   28-Mar-2022 :   Fixed double correcting the angle.
#   04-Apr-2022 :   Wrapped code into class
#   11-Apr-2022 :   Greatly improved documentation on Getstars and added
#                   some other functions.
#   17-Apr-2022 :   Implementation of finding fields/subfields locations
#   18-Apr-2022 :   Added __init__ documentation
#   20-Apr-2022 :   Fixed bug where lists were linking, first successful
#                   run. Renamed file from ukirtSurot.py to createCMD.py
#   21-Apr-2022 :   Starts reddening vector calculation
#   22-Apr-2022 :   Finished reddening vector function.
#-----------------------------------------------------------------------

#Importing all necessary packages
import os
import pickle
from astropy import coordinates as coord
from astropy import units as u
from matplotlib.pyplot import hist
import numpy as np
import parameters as pram
import ipdb

class cmd:
    def __init__(self, ra=None, dec=None, l=None, b=None, year = pram.year, edge_length=0.25, findvec = False, fieldType = 'field', field_ind=0, rc_dict={}):
        """The cmd class creates cmds at a given location. The class can be fed the location of a square area, denoted 
        as pixel, to find the cmd of or use the field and subfields themselves. The reddening vector of a given field
        or subfield is able to be calculated as well. If the reddening vector isn't being found then ra and dec
        coordinates are necessary to be given. Additionally, the galactic coordinates l and b can be either given
        to class or calculated internally based on the ra and dec.

        Parameters
        ----------
        ra : float, optional
            The right ascension of the pixel center, by default None
        dec : float, optional
            The declination of the pixel center, by default None
        l : float, optional
            The galactic longitude of the pixel center, by default None
        b : float, optional
            The galactic latitude of the pixel center, by default None
        year : int, optional
            The year of desired observation period, by default pram.year
        edge_length : float, optional
            The edge length of a given pixel, by default 0.25
        findvec : bool, optional
            Whether the reddening vector is to be calculated, by default False
        fieldType : str, optional
            When finding the reddening vector, this denotes the size of
            area to use. Can be 'field' or 'subfield', by default 'field'
        field_ind : int, optional
            The index of the field or subfield that is to be used. Values can be found in
            the pdf files 'ukirt_fieldgrid_2017.pdf' or 'ukirt_subfieldgrid_2017.pdf', by default 0
        """

        # Here we start declaring variables as global

        self.year = year
        self.edge_length = edge_length

        self.findvec = findvec
        self.fieldType = fieldType
        
        # Obtaining the UKIRT field positions
        self.fieldData = self.readUKIRTfields()

        # Checks if we are finding the reddening vector
        # This will change out we find our stars if so and
        # which limiting dictionaries we shall use.
        # There are two limiting dictionaries to be used. One for 
        # stellar population limits and the other for color-magnitude cuts.
        if self.findvec == True:

            # Retrieve list of Field/Subfield locations
            fieldRange = self.fieldLocations()
            
            # Extract ra/dec min/max values for specific field denoted by
            # field_ind where the integer corresponds with the value of field
            # which can be found in the figures 'ukirt_subfieldgrid_2017.pdf'
            # or 'ukirt_fieldgrid_2017.pdf' in the 'figs/' folder by default.
            self.ra_min = fieldRange[field_ind][1]
            self.ra_max = fieldRange[field_ind][2]
            self.dec_min = fieldRange[field_ind][3]
            self.dec_max = fieldRange[field_ind][4]

            # Since when choosing fields we assume ra and dec isn't given, calculate it
            self.ra = (self.ra_max + self.ra_min)/2
            self.dec = (self.dec_max + self.dec_min)/2

            #Checks if l and b are given. If not, calculate them
            if type(l)==type(None) or type(b)==type(None):
                c = coord.SkyCoord(ra=self.ra*u.degree,dec=self.dec*u.degree ,frame='icrs')
                l = c.galactic.l.degree
                b = c.galactic.b.degree
        
            self.l = l
            self.b = b
            
            self.raw_field_inds = self.findFields(plot_fields=False)
            self.limit_dict = {'altMAD':[-0.1,0.1], 'MAD':[-0.1,0.1]}
            #self.limit_dict = {'offset':[-0.1,0.1]}
            #self.limit_dict={}
            self.cm_dict = {}
            
        else:

            # Set our ra and dec
            self.ra = ra
            self.dec = dec

            #Checks if l and b are given. If not, calculate them
            if type(l)==type(None) or type(b)==type(None):
                c = coord.SkyCoord(ra=self.ra*u.degree,dec=self.dec*u.degree ,frame='icrs')
                l = c.galactic.l.degree
                b = c.galactic.b.degree
        
            self.l = l
            self.b = b

            self.ra_max = self.ra+self.edge_length/2.
            self.ra_min = self.ra-self.edge_length/2.
            self.dec_max = self.dec+self.edge_length/2.
            self.dec_min = self.dec-self.edge_length/2.

            self.raw_field_inds = self.findFields(plot_fields=False)

            self.limit_dict = {'N':[10,1e6],'altN':[3,1e6],'offset':[-0.1,0.1]}
            self.cm_dict = {'altmag':[12,16],'delta':[-1,5]}

        
        self.getStars(limit_dict=self.limit_dict)
        self.color_mag_cut(self.cm_dict,percentile_cut=True)

        self.coeffs = self.calcReddeningVec(rc_dict)
        self.red_vec = self.coeffs[0]


    def readUKIRTfields(self, filename='ukirtFieldLocations.sav', dir='fieldLocations/'):
        """Reads in the field data for the United Kindgon InfraRed Telescope to give
        a list of dictionaries containing the field data.

        Parameters
        ----------
        filename : str, optional
            Name of the .sav file hosting the data, by default 'ukirtFieldLocations.sav'
        dir : str, optional
            Directory where the file is located, by default 'fieldLocations/'

        Returns
        -------
        fieldData : list
            List of dictionaries which contains each subfield and all its corresponding information.
        """
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
        """Function which gives the RA and DEC minimum and maximum values for a specific ccd.
        See fieldLocations/ukirtFieldLocations.txt for labeling of the fields.

        Parameters
        ----------
        field : string
            The label for interested field.
        ccd : string
            Desired ccd, labeled 1 through 4.
        year : string, optional
            String for the year of desired observation period, by default pram.year

        Returns
        -------
        _type_
            _description_
        """
        # Each ccd has fov of 13.6'x13.6', so each field has 54.4'x54.4'
        
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


    def fieldLocations(self):
        """
        Generates the minimum and maximum RA and DEC values for the fields or
        subfields when finding the reddening vector.

        Returns
        -------
        fieldRange : list
            List containing the index of field including min and max values.
        """
        
        vecfieldData = []

        # Creates and appends new list to hold only data in our observation year
        newfieldData = []
        for i in range(len(self.fieldData)):
            if int(self.fieldData[i]["year"]) == int(self.year):
                newfieldData.append(self.fieldData[i].copy())

        if self.fieldType == 'field':
            # This takes the field generated by all four ccd's in all 4 positions
            vecfieldData = [newfieldData[i:i+16] for i in range(0, len(newfieldData), 16)]

        elif self.fieldType == 'subfield':
            # Takes the subfield made by combining one location from each ccd
            tempfield = [newfieldData[i:i+16] for i in range(0, len(newfieldData), 16)]
            for field in tempfield:
                for i in range(0,4,1):
                    # Grabs every 4th field
                    vecfieldData.append(list(np.array(field)[[i,i+4,i+8,i+12]]))

        else:
            exit("ERROR: Invalid field type or no field type selected.")
        
        # Restructuring
        ordered_fields = []

        for dict_list in vecfieldData:

            # Creates a dictionary with the framework of the first entry
            all_dict=dict_list[0]

            for i in np.arange(len(dict_list)-1)+1:
                # Loops over all keys
                for key in all_dict.keys():
                    # Appends the dictionary key with the value of the next datapoint of the i-th star
                    element = np.append(all_dict[key],dict_list[i][key])
                    all_dict[key] = element
            
            ordered_fields.append(all_dict)

        # List which will hold field number and min/max values
        fieldRange = []
        
        for i in range(len(ordered_fields)):

            # Calculating the min and max values 
            ramin = np.amin(ordered_fields[i]['RAmin'])
            ramax = np.amax(ordered_fields[i]['RAmax'])
            decmin = np.amin(ordered_fields[i]['DECmin'])
            decmax = np.amax(ordered_fields[i]['DECmax'])

            # Putting into the list
            fieldRange.append([i,ramin,ramax,decmin,decmax])

        return fieldRange


    def findFields(self,new_ra=None,new_dec=None,plot_fields=True):
        """Gives the indexes of fields from the fieldData list of UKIRT fields for which
        our CMD will use. It does so by checking the range of coordinates for each field and
        seeing if our desired location lands within it.

        Parameters
        ----------
        new_ra : float, optional
            New ra value to update with, by default None
        new_dec : float, optional
            New dec value to update with, by default None

        Returns
        -------
        raw_field_inds : list
            List of indices which the pixel lands in.
        """

        # Checks if new RA and DEC coords were given
        if new_ra == None and new_dec == None:
            ra_max = self.ra_max
            ra_min = self.ra_min
            dec_max = self.dec_max
            dec_min = self.dec_min
        else:
            # Calculates max and min if new coords given
            ra_max = new_ra + self.edge_length/2.
            ra_min = new_ra - self.edge_length/2.
            dec_max = new_dec + self.edge_length/2.
            dec_min = new_dec - self.edge_length/2.
        
        # Creates empy list
        raw_field_inds = []
        for i in range(len(self.fieldData)):

            # Make sure we use the year we want
            if self.fieldData[i]['year']==str(self.year):
                
                # Checks if the field lies within the RA range
                if (ra_min<=self.fieldData[i]['RAmax'] and ra_min>=self.fieldData[i]['RAmin']) or \
                        (ra_max<=self.fieldData[i]['RAmax'] and ra_max>=self.fieldData[i]['RAmin']) or\
                        (ra_max>=self.fieldData[i]['RAmax'] and ra_min<=self.fieldData[i]['RAmin']):

                    # Checks if the field lies in the DEC range after confirming it's in RA range.
                    if (dec_min<=self.fieldData[i]['DECmax'] and dec_min>=self.fieldData[i]['DECmin']) or \
                            (dec_max<=self.fieldData[i]['DECmax'] and dec_max>=self.fieldData[i]['DECmin']) or \
                            (dec_max>=self.fieldData[i]['DECmax'] and dec_min<=self.fieldData[i]['DECmin']):
                    
                        # Append the fieldData index i to our list
                        raw_field_inds.append(i)
                        
        return raw_field_inds


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
        self.allStarDict : dict
            Dictionary of all stars within the field range, no filters applied.
        self.filterStarDict : dict
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

        
        self.allStarDict = all_dict
        self.filterStarDict = filt_dict


    def color_mag_cut(self,limit_dict=None,percentile_cut=False):
        """Function which performs a color-magnitude cut on the
        filtered star dictionary. Creates a cmd class dictionary
        which contains the stars passing the cut, called self.fitStarDict.

        Parameters
        ----------
        limit_dict : dict, optional
            A dictionary containing filters on any key found within the stellar
            data. Keys of limit_dict must match keys found in the data files
            and provide lower and upper bound ranges. Primary filters include on 
            number of stars or the magnitude of stars.
            Ex: limit_dict = {'altmag':[12.16],'N':[10,1e6]}
            Defaults to None
        percentile_cut : bool, optional
            Whether to perform a percentile cut in (H-K) color
            where only stars above 15th percentile are kept, by default False
        """

        # Checks if there is a limit dictionary set
        if type(limit_dict)==type(None):
            print('No limits on color or magnitude')

        # Make a boolean index array to initialize our color_mag cuts
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


    def updateCMD(self,new_ra,new_dec,l=None,b=None):
        """Updates the cmd if location changes in calculations.

        Parameters
        ----------
        new_ra : float, optional
            New ra value to update with, by default None
        new_dec : float, optional
            New dec value to update with, by default None
        l : float, optional
            The galactic longitude of the pixel center, by default None
        b : float, optional
            The galactic latitude of the pixel center, by default None
        """

        # Set fundamental box properties, center and length of box
        self.ra = new_ra
        self.dec = new_dec
        if type(l)==type(None) or type(b)==type(None):
            c = coord.SkyCoord(ra=new_ra*u.degree,dec=new_dec*u.degree ,frame='icrs')
            l = c.galactic.l.degree
            b = c.galactic.b.degree
        self.l = l
        self.b = b
        
        # Set the upper and lower limits for our box in ra and dec
        self.ra_max = self.ra+self.edge_length/2.
        self.ra_min = self.ra-self.edge_length/2.
        self.dec_max = self.dec+self.edge_length/2.
        self.dec_min = self.dec-self.edge_length/2.
        
        # Gathers fields and gets stars within them, performing cuts
        self.raw_field_inds = self.findFields()
        self.getStars(skip_read=True,limit_dict=self.limit_dict)
        self.color_mag_cut(self.cm_dict,percentile_cut=True)
        self.histMag(plot=False)
        
        print('CMD updated!')


    def calcReddeningVec(self, limit_dict=None):


        # Creating a dictionary which will house our suspected RC stars
        self.RC_dict = {}
        
        # Starting off with our good indices being those within coord range
        good_inds = np.ones(len(self.filterStarDict['RA']),dtype=bool)

        # Looping over keys
        for key in self.filterStarDict.keys():

            # Checks if the specified key is within our limiting filter dictionary.
            if key in limit_dict.keys():

                # Gathers indices where the stars in the dictionary pass the filter.
                low_inds = self.filterStarDict[key]>=limit_dict[key][0]
                high_inds = self.filterStarDict[key]<=limit_dict[key][1]

                # Gathers all elements that don't have "nan" entries
                try:
                    # Returns False if 'nan' and True if not, basically opposite of Numpy.isnan
                    nan_inds = not (np.isnan(self.filterStarDict[key]))
                except:
                    # If there are no 'nan' values, return a list of all True
                    nan_inds = np.ones(len(self.filterStarDict[key]),dtype=bool)
                    pass

                # Uses logical 'and' to gather all indices that have passed checks
                good_inds = good_inds & low_inds & high_inds & nan_inds

            # If no filter key, just checks for 'nan' same way as above
            else:
                try:
                    good_inds = good_inds & (not np.isnan(self.filterStarDict[key]))
                except:
                    nan_inds = np.ones(len(self.filterStarDict[key]),dtype=bool)
                    good_inds = good_inds & nan_inds

        # Appends our filtered dictionary with values in 'good_inds'
        for key in self.filterStarDict.keys():
            self.RC_dict[key] = self.filterStarDict[key][good_inds]

        # Here we take the density of points and find where it is max
        # This is done by splicing the 2-D space into a bunch of slices
        # and creating a histogram for each slice, finding the peak of each.
        # The (H-K) color and mag K is saved for each peak to a list of points
        # to be used for the fit.
        ipdb.set_trace()
        # Splitting the boxed area in cmd-space to make histograms in each color line
        # Need to add 1 extra value so that last one can be cut off, since we loop using
        # Starting location + width
        histLocations = np.linspace(limit_dict['delta'][0], limit_dict['delta'][1], 51)
        # Width of bins
        width = round(abs(histLocations[0]-histLocations[1]), 6)
        # Removing the last point, giving us 50 bins!
        histLocations = histLocations[:-1]
        
        # Empty dictionary for the fit data
        vecFitList = {'delta':[],'altmag':[]}
        
        # Loops over each of these color lines
        for xbin in histLocations:
            # Resets the good indices
            good_inds = np.ones(len(self.RC_dict['RA']),dtype=bool)

            # temporary dictionary
            tempDict = {'delta':[],'altmag':[]}

            # Sees if the values are within the bin range and limits indices to those
            low_inds = self.RC_dict['delta']>=xbin
            high_inds = self.RC_dict['delta']<=(xbin+width)
            good_inds = good_inds & low_inds & high_inds

            # Adds points to the temp dict for each key
            for key in tempDict.keys():
                tempDict[key] = self.RC_dict[key][good_inds]

            # Calculates the histogram
            hist, bin_edges = np.histogram(tempDict['altmag'], bins=50)

            # Gets index of maximum histogram value and gets value of magnitude at max
            index = np.argmax(hist)
            maxDensity = (bin_edges[index]+bin_edges[index+1])/2
            # Adds color and magnitude to fit list
            vecFitList['delta'].append(xbin)
            vecFitList['altmag'].append(maxDensity)
            
        # Fits a line to the max density points
        coeffs = np.polyfit(vecFitList['delta'],vecFitList['altmag'],1)
        self.vecFitList = vecFitList
        # Returns the coefficients
        return coeffs

if __name__ == '__main__':
    rc_dict={'altmag':[12, 16], 'delta':[1.0, 2.0], 'altMAD':[-0.1,0.1], 'MAD':[-0.1,0.1]}
    cmd_test = cmd(findvec=True, fieldType='subfield', field_ind=24, rc_dict=rc_dict)
    print(cmd_test.coeffs[0])
