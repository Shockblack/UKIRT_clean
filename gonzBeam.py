#-------------------------------------------------------
# File name: gonzBeam.py
# 
# Creates a data list for Gonzalez data from
# a BEAM dat file from Gonzalez et al. 2012.
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
#-------------------------------------------------------

# Importing necessary packages and files
import numpy as np #cause you know, numpy
import ipdb #For error testing

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


    def writeFile(self):

        self.gonzData = []

        # Creates one large list for data organization
        for i in range(len(self.rawGonzDat)):
            data = self.coords[i] + self.gonzAK[i] + self.rawGonzDat[i][2:]
            self.gonzData.append(data)

        #Opens file path
        fileout = open(self.path+'gonzData.dat', 'w')

        #Writes data to file with organization of
        # ['ra', 'dec', 'l', 'b', 'AK', 'sigmaAK', 'E_JK', 'sigmaE_JK', '[Fe/H]']
        fileout.write('# ra, dec, l, b, AK, sigmaAK, E_JK, sigmaE_JK, [Fe/H] \n')
        for i in range(len(self.rawGonzDat)):

            for d in self.gonzData[i]:
                fileout.write(str(d) +' ')
            
            fileout.write('\n')
        
        fileout.close()




if __name__=='__main__':

    ipdb.set_trace()
    gonz = gonzDat() #Initiates class. Edit 'path' and 'filename' parameters to change path
    gonz.calcRaDec()
    gonz.calcAK()
    gonz.writeFile()
