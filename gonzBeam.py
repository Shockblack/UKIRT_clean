#-------------------------------------------------------
# File name: gonzBeam.py
# 
# Creates a data list for Gonzalez data from
# a BEAM dat file from Gonzalez et al. 2012 and 2018.
# File contains l, b, E(J-K), sig_E(J-K), and
# [Fe/H]. We primarily use the first four entries
# for our work.
#
# Programmer: Aiden Zelakiewicz (zelakiewicz.1@osu.edu)
#
# For questions, reach out to:
#   Aiden Zelakiewicz   (zelakiewicz.1@osu.edu)
#   Samson Johnson      (johnson.7080@osu.edu)
#
# Revision History:
#   06-Feb-2022: File Created
#   14-Feb-2022: Added Ra and Dec calculations, moved
#                file writing to its own function.
#                Changed delimiter in file to ',' rather
#                than spaces.
#   19-Mar-2022: Added function to make the Gonzalez data
#                file coordinates line up with UKIRT
#                coords and export to different file.
#-------------------------------------------------------

# Importing necessary packages and files
import numpy as np #cause you know, numpy
import ipdb #For error testing
import pandas as pd #Used to create comparison file

# Importing astropy to calculate ra and dec
# This is to make gonz data file a bit easier to use
# with existing map data
from astropy import coordinates as coord
from astropy import units as u

class gonzDat:

    def __init__(self,path='../data/gonz/',filename='EJK_beam.dat'):
        
        self.path = path
        self.filename = filename

        #Data has organization of ['l','b','E_JK','sigmaE_JK','[Fe/H]']
        self.rawGonzDat = np.genfromtxt(self.path+self.filename).tolist()


    def calcRaDec(self, edgelength=0.025):

        self.edgelength = edgelength # 0.025 by 0.025 deg^2 is the pixel size of Gonz. 2012

        self.coords = []

        # probably a better way to do this, empty l b lists
        lList = []
        bList = []
        
        # Creates a location list of l and b values
        for i in range(len(self.rawGonzDat)):
            lList.append(self.rawGonzDat[i][0])
            bList.append(self.rawGonzDat[i][1])
        
        # creates SkyCoord object taking in l and b in
        c = coord.SkyCoord(l=lList*u.degree,b=bList*u.degree ,frame='galactic')

        # Calculating Ra and Dec
        raList = c.icrs.ra.degree
        decList = c.icrs.dec.degree

        for i in range(len(lList)):
            coordList = [raList[i],decList[i],lList[i],bList[i]]
            self.coords.append(coordList)


    def calcAK(self):

        self.gonzAK = []

        for i in range(len(self.rawGonzDat)):

            # Using conversion values Gonzalez et al. 2012 used from Cardelli et al. 1989
            # AK = 0.689 * E(J-K)
            ak = 0.689*self.rawGonzDat[i][2]
            akerr = 0.689*self.rawGonzDat[i][3]

            data = [ak,akerr]

            self.gonzAK.append(data)


    def writeFile(self, filename='gonzData.dat'):

        self.gonzData = []

        # Creates one large list for data organization
        for i in range(len(self.rawGonzDat)):
            data = self.coords[i] + self.gonzAK[i] + self.rawGonzDat[i][2:]
            self.gonzData.append(data)

        #Opens file path
        fileout = open(self.path+filename, 'w')

        #Writes data to file with organization of
        # ['ra', 'dec', 'l', 'b', 'AK', 'sigmaAK', 'E_JK', 'sigmaE_JK', '[Fe/H]']
        fileout.write('# ra, dec, l, b, AK, sigmaAK, E_JK, sigmaE_JK, [Fe/H] \n')
        for i in range(len(self.rawGonzDat)):

            for d in range(len(self.gonzData[i])):
                if d != 8:
                    fileout.write(str(self.gonzData[i][d]) + ',')
                else:
                    fileout.write(str(self.gonzData[i][d]))
            
            fileout.write('\n')
        
        fileout.close()

    def compare_map(self, map1='../data/gonz/gonzData.dat', map2='maps/map_PSF_2017_2_gonzGrid.map', \
                    l_axis = 2, b_axis = 3, fileout='gonzDataScaled.dat'):
        """
        Takes the data from two seperate maps and compares
        their galactic locations to scale the map files to be of same length. 
        Saves this new data as a seperate csv file. Check function to
        make sure pd.read_csv has the correct kwargs for files using.

        Code inspirations and techniques used:
        https://moonbooks.org/Articles/How-to-delete-rows-with-values-below-and-above-a-minimum-and-maximum-value-in-a-pandas-data-frame-/
        https://stackoverflow.com/questions/62519791/finding-duplicates-in-two-dataframes-and-removing-the-duplicates-from-one-datafr
        https://stackoverflow.com/questions/17978133/python-pandas-merge-only-certain-columns

        Dependencies
        ------------
        pandas

        Parameters
        ----------
        map1 : string
            Filepath of 'larger' map data to be edited.
        map2 : string
            Filepath of map to have locations matched to.
        l_axis : int
            Axis on map1 where l data is located.
        b_axis : int
            Axis on map1 where b data is located.
        fileout : string
            Name of file resulting data to be outputted to.
            Goes to default location of self.path if a seperate
            path is not specified.
        """

        #Loads both the map files
        comp_map = pd.read_csv(map2,header=None)
        large_map = pd.read_csv(map1,skiprows=1, header=None)

        # Only keeps values of larger map which have the same location as the smaller/comparison map
        large_map = pd.merge(large_map, comp_map[[l_axis,b_axis]], on=[l_axis,b_axis], how='inner')

        # Organizes the map by decreasing l and b values
        large_map = large_map.sort_values(by=[l_axis,b_axis], ascending=False)
        
        # saves map to csv file, doesn't include header
        large_map.to_csv(self.path+fileout, header=False, index=False)


if __name__=='__main__':
    # EJK_BEAM.dat is 2012 data and EJK_BEAM_2018.txt is 2018 gonz data
    gonz = gonzDat(filename='EJK_BEAM_2018.txt') #Initiates class. Edit 'path' and 'filename' parameters to change path
    gonz.calcRaDec()
    gonz.calcAK()
    ipdb.set_trace()
    gonz.writeFile(filename='gonzData_2018.dat')
    gonz.compare_map(map1='../data/gonz/gonzData_2018.dat', fileout='gonzDataScaled_2018.dat')