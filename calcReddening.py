#-----------------------------------------------------------------------
# File name: calcReddening.py
# 
# Calculates the reddening and extinction using the method by Surot et al.
# 2020 and derives code from predecessor file "mapper.py." Meant to do so
# by taking a given field/subfield and dividing it up into smaller "pixels."
# this will eliminate the need to worry about points overlapping fields where
# derived reddening vectors are different.
#
#-----------------------------------------------------------------------

# Importing other Python files from project
import createCMD

# Importing necessary packages
import numpy as np
import parameters as pram
from astropy import coordinates as coord
from astropy import units as u
import ipdb

class mapper:
    def  __init__(self, cmdclass, year=pram.year, field_ind=0, num_pixels=128):

        # Setting the variables as global variables
        self.year = year
        self.field_ind = field_ind

        # This is the class from createCMD.py which stores values such as max/min ra/dec
        # This also includes a dictionary of the stars in the RC, named RC_dict
        self.cmd = cmdclass

        self.raw_fieldData = self.cmd.readUKIRTfields()
        
        # Select by year
        self.fieldData = []
        for field in self.raw_fieldData:
            if int(field['year']) == year:
                self.fieldData.append(field)

        # Now we must find the size/location of the fields
        # by taking the field max/min and dividing by num of pixels
        self.edge_length = abs(self.cmd.ra_max - self.cmd.ra_min) / num_pixels

        self.dimension = self.findDimension()

    def gen_grid(self, grid_grow_scale=0.025):
        
        ra_span = self.cmd.ra_max-self.cmd.ra_min
        dec_span = self.cmd.dec_max-self.cmd.dec_min

        ra_max = self.cmd.ra_max+grid_grow_scale*ra_span
        ra_min = self.cmd.ra_min-grid_grow_scale*ra_span

        dec_max = self.cmd.dec_max+grid_grow_scale*dec_span
        dec_min = self.cmd.dec_min-grid_grow_scale*dec_span

        ra_range = np.arange(ra_min, ra_max, self.edge_length)
        dec_range = np.arange(dec_min, dec_max, self.edge_length)


        self.grid_pixel_centers = []

        #Calculates RA and DEC right when map files begins being made for earier use and to save time.
        for ra in ra_range:
            for dec in dec_range:
                c = coord.SkyCoord(ra=ra*u.degree,dec=dec*u.degree ,frame='icrs')
                l = c.galactic.l.degree
                b = c.galactic.b.degree
                self.grid_pixel_centers.append([ra,dec,l,b])

        return


    def getNABStars(self, pixel_limits):

        # Creating empy dictionary
        pixel_RC_dict = {}
        
        # Getting all the stars within a specific pixel
        good_ra_high = self.cmd.RC_dict['RA']<=pixel_limits['RAmax'] 
        good_ra_low = self.cmd.RC_dict['RA']>=pixel_limits['RAmin'] 
        good_ra_inds = good_ra_low & good_ra_high

        good_dec_high = self.cmd.RC_dict['Dec']<=pixel_limits['DECmax']
        good_dec_low = self.cmd.RC_dict['Dec']>=pixel_limits['DECmin'] 
        good_dec_inds = good_dec_low & good_dec_high 
        
        good_inds = good_ra_inds & good_dec_inds

        for key in self.cmd.RC_dict.keys():
            pixel_RC_dict[key] = self.cmd.RC_dict[key][good_inds]

        # First going to restructure and reorganize because
        # I can only think in very specific ways (sorry)
        # While doing this also going to calculate the distance

        # Empty lists
        RC_dict_org = []
        dist_list = []

        for i in range(len(pixel_RC_dict['RA'])):
            # Temporary dictionary
            temp_dict = {}
            # Looping over the keys to make the headers
            for key in pixel_RC_dict.keys():
                temp_dict[key] = pixel_RC_dict[key][i]
            # Just creating some variables for easier use
            ra = self.cmd.ra
            dec = self.cmd.dec
            # Caclulating the distance
            distance = np.sqrt((temp_dict['RA']-ra)**2 + (temp_dict['Dec']-dec)**2)
            temp_dict['dist'] = distance
            dist_list.append(distance)
            RC_dict_org.append(temp_dict)
        
        dist_array = np.array(dist_list)

        # Find the indices of the 20 clostest stars
        ind_array = np.argsort(dist_array)

        # Saving the max distance of these 20 stars in 
        # order to use in calculating weights later on
        self.max_dist = max(dist_array[ind_array][ : 20])

        # Getting the 20 stars
        self.NABStars = np.array(RC_dict_org)[ind_array][ : 20]

        # Creating a list of colors for finding Median Absolute Deviation (MAD)
        color_array = []
        for star in self.NABStars:
            color_array.append(star['delta'])
        
        #color_array = np.array(color_array)
        
        # Now we calculate the MAD and remove any
        # stars with color 3 MAD away from the median
        median_color = np.median(color_array)
        abs_difference = np.absolute(np.array(color_array) - median_color)
        self.mad = np.median(abs_difference)
        
        # Finding indices of stars 3 MAD away
        bad_ind = []
        # Looping over all stars
        for i in range(len(color_array)):
            # Checking if 3 MAD away
            if np.abs(color_array[i]-median_color) > 3*self.mad:
                # Appending bad index
                bad_ind.append(i)
        
        # Deleting the stars over 3 MAD away
        self.NABStars = np.delete(self.NABStars, bad_ind)


    def calcWeights(self):

        # Empty list to store the weight values
        self.weights = []

        # Caclculating the weights for each star
        for star in self.NABStars:

            # The ratio of the distance to the largest one
            dist_ratio = star['dist']/self.max_dist

            # The ith weight
            star['weight'] = 1.0 - 0.75 * dist_ratio**2

            # Appending to the list of weights
            self.weights.append(star['weight'])
        
        # Calculating the sum in order to normalize
        # the weights to sum to 1
        weight_sum = np.sum(self.weights)

        # Normalizing the weight list
        self.weights = np.array(self.weights)/weight_sum

        # Going to loop over again to assign the normalized
        # weights to the dictionaries
        for star in self.NABStars:

            star['weight'] = star['weight']/weight_sum


    def calcReddening(self):

        # Creating an empty variable
        average_color = 0.

        # Summing the color excess with their weights
        for star in self.NABStars:
            average_color += star['weight'] * star['delta']

        # Dividing by number of points to get the average
        color = average_color / np.sum(self.weights)

        # Calculating color E(H-K) from E(H-K) = (H-K) - (H-K)_0
        # where we get (H-K)_0 from Nataf et al. 2021 and is similar
        # to the value of our data, so we decide to use Nataf's value.
        reddening = color - 0.15

        # Also calculate the extinction using the reddening vector
        # found in the file createCMD.py
        extinction = self.cmd.red_vec * reddening
        
        return reddening, extinction

    def findDimension(self, num_bins=10):
        """Function which will adjust the number of pixels within our map until 
        there is a minimum of 20 stars within every pixel. It will raise or lower
        the amount of points until this minimum is hit. This takes a dimension, n, 
        and splits the RA space into n "lines." In each of these lines, a 1D histogram
        is applied to it using Numpy. The minimum is recorded and compared against the
        other lines until an absolute minimum is acquired. If it is above the threshold
        of 20 stars, the dimension is increased by 1 and ran over again. This happens
        until the minimum drops below 20 and the dimension is recorded.

        Parameters
        ----------
        num_bins : int, optional
            The starting amount of pixels to begin looping over, by default 50

        Returns
        -------
        dimension : int
            The final dimension the subfield is broken up into. The amount of pixels
            used is therefore dimension * dimension.
        """
        
        # Set starting minimum star value so that is is global
        min_stars = 10000
        # Adding a counter to know how many iterations it took
        counter = 0

        # Looping as long as the amount of stars is 20 or greater
        while min_stars >= 20:

            # Splitting the subfield's RA space into lines (bins)
            # Need 1 extra so that we can clip last one and not go over box edge
            histLocations = np.linspace(self.cmd.ra_min, self.cmd.ra_max, num_bins+1)
            # Width of bins
            width = abs(histLocations[0]-histLocations[1])
            # Clipping last point
            histLocations = histLocations[:-1]

            # Loops over each of these RA lines
            for xbin in histLocations:
                
                # Resets the good indices
                good_inds = np.ones(len(self.cmd.RC_dict['RA']),dtype=bool)

                # temporary dictionary
                tempDict = {'RA':[],'Dec':[]}

                # Sees if the values are within the bin range and limits indices to those
                low_inds = self.cmd.RC_dict['RA']>=xbin
                high_inds = self.cmd.RC_dict['RA']<=(xbin+width) # This width is why we have the extra bin and clip it
                good_inds = good_inds & low_inds & high_inds

                # Adds points to the temp dict for each key
                for key in tempDict.keys():
                    tempDict[key] = self.cmd.RC_dict[key][good_inds]

                # Calculates the histogram
                hist, bin_edges = np.histogram(tempDict['Dec'], bins=num_bins)

                # Records the local minimum
                line_min = np.amin(hist)
                # Checks if it's a local min or global min
                if line_min < min_stars:
                    # If global min, set it to the variable
                    min_stars = line_min
            
            print("Min Stars =",min_stars)
            print("Num bins = ",num_bins)
            # Increase bin count and counter
            num_bins += 1
            counter += 1
        
        # Need to remove 2 bins, since if we are less than 20, we need to remove 1 from the 
        # addition at the end and another to go to the dimension before we dropped below.
        dimension = num_bins - 2
        print("Program converged on ",dimension,"x",dimension," pixels after ",counter," iterations.")
        
        return dimension

    def pixelInfo(self, dimension=None):
        
        # Checking if a custom dimension isn't given
        # If none is given, uses the automatic one 
        # calculated from function findDimension()
        if type(dimension) == type(None):
            dimension = self.dimension
        
        # Creating the dictionary to house the pixel info
        self.pixel_data = {'RA' : [], 'Dec' : [], 'l' : [], 'b' : [], 'edgelength' : [], 'RAmax' : [], \
            'RAmin' : [], 'DECmax' : [], 'DECmin' : [], 'reddening' : [], 'extinction' : [], }

        # Splitting field into the pixels
        
        # Doing RA
        rabins = np.linspace(self.cmd.ra_min, self.cmd.ra_max, dimension+1)
        rabins = rabins[:-1]

        # and Dec
        decbins = np.linspace(self.cmd.dec_min, self.cmd.dec_max, dimension+1)
        decbins = decbins[:-1]

        # Start looping over RA and Dec space
        # Begin with dec first
        for i in range(len(decbins)):
            # Assign one edge of bin
            dec_min = decbins[i]
            # Set other edge to either next bin or the max value of field if last point
            if i == dimension-1:
                dec_max = self.cmd.dec_max
            else:
                dec_max = decbins[i+1]
            
            # Start looping over ra
            for k in range(len(rabins)):
                # Using one edge of bin
                ra_min = rabins[k]
                # Same thing, checking for edge case
                if k == dimension-1:
                    ra_max = self.cmd.ra_max
                else:
                    ra_max = rabins[k+1]
                
                # Calculating the ra/dec and l/b of pixel center
                # ra and dec
                ra = (ra_min+ra_max)/2
                dec = (dec_min+dec_max)/2

                # l and b
                c = coord.SkyCoord(ra=ra*u.degree,dec=dec*u.degree ,frame='icrs')
                l = c.galactic.l.degree
                b = c.galactic.b.degree

                # Giving the dictionary location values so stars can be grabbed
                self.pixel_data['RA'].append(ra)
                self.pixel_data['Dec'].append(dec)
                self.pixel_data['l'].append(l)
                self.pixel_data['b'].append(b)
                self.pixel_data['edgelength'].append(np.abs(ra_max-ra_min))
                self.pixel_data['RAmin'].append(ra_max)
                self.pixel_data['RAmax'].append(ra_min)
                self.pixel_data['DECmin'].append(dec_max)
                self.pixel_data['DECmax'].append(dec_min)

                # Gathering a dict of pixel location limits for getNABStars
                pixel_limits = {'RAmax' : ra_max, 'RAmin' : ra_min,\
                    'DECmax' : dec_max, 'DECmin' : dec_min}

                # Get all the stars within a pixel
                self.getNABStars(pixel_limits)

                # Calculate the reddening and extinction for the specific pixel
                self.calcWeights()
                reddening, extinction = self.calcReddening()

                # Appending reddening and extinction values to dict
                self.pixel_data['reddening'].append(reddening)
                self.pixel_data['extinction'].append(extinction)

        # Calculating what the edgelength of the map is
        self.edge_length = self.pixel_data['l'][0] - self.pixel_data['l'][1]
        

    def saveMap(self, path = 'maps/'):
        header = 'RA,DEC,l,b,edgelength,N_stars,A,Aerr,B,Berr,M_RC,M_RCerr,sigma_RC,sigma_RCerr,N_RC,N_RCerr,FinalColorRC,RealSTDcolorOpt,STDtotal,VarianceMin,MU1opt,MU2opt'
        #date 
        #year edge_length number_pixels 
        #header for columns
        #we will probably want some header info
        with open(path+'field_pixel_data_'+str(self.field_ind)+'.map','w') as outfile:

            for i in range(len(self.pixel_data['RA'])):
                write_string = ""

                for key in self.pixel_data.keys():
                    write_string+='%0.6f,'%self.pixel_data[key][i]

                outfile.write(write_string[:-1]+'\n')
        return


if __name__ == "__main__":
    import pandas as pd

    RC_limits_df = pd.read_csv('RC_limits_byEye.csv')
    
    for i in range(56):
        field_limits = RC_limits_df.iloc[i]
        field_ind = i

        rc_dict={'altmag':[field_limits['Kmin'], field_limits['Kmax']],\
            'delta':[field_limits['HKmin'], field_limits['HKmax']], 'altMAD':[-0.1,0.1], 'MAD':[-0.1,0.1]}
        cmd_test = createCMD.cmd(findvec=True, fieldType='subfield', field_ind=field_ind, rc_dict=rc_dict)
        map = mapper(cmd_test, field_ind=field_ind)
        map.pixelInfo()
        map.saveMap()