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
#-------------------------------------------------------

# Importing necessary packages and files
import numpy as np #cause you know, numpy
import ipdb #For error testing

class gonzDat:

    def __init__(self,path='../data/gonz/',filename='EJK_beam.dat'):
        
        self.path = path
        self.filename = filename

        #Data has organization of ['l','b','E_JK','sigmaE_JK','[Fe/H]']
        self.rawGonzDat = np.genfromtxt(self.path+self.filename).tolist()

    def calcAK(self):

        self.gonzData = []

        fileout = open(self.path+'gonzData.dat', 'w')

        for i in range(len(self.rawGonzDat)):

            # Using conversion values Gonzalez et al. 2012 used from Cardelli et al. 1989
            # AK = 0.689 * E(J-K)
            ak = 0.689*self.rawGonzDat[i][2]
            akerr = 0.689*self.rawGonzDat[i][3]

            # Appending the extinction values to give 
            # ['l', 'b', 'E_JK', 'sigmaE_JK', '[Fe/H]', 'AK', 'sigmaAK']

            data = self.rawGonzDat[i]
            data.append(ak)
            data.append(akerr)

            self.gonzData.append(data)

            for d in self.gonzData[i]:
                fileout.write(str(d) +' ')
            
            fileout.write('\n')
        
        fileout.close()



if __name__=='__main__':

    ipdb.set_trace()
    gonz = gonzDat() #Initiates class. Edit 'path' and 'filename' parameters to change path
    dat = gonz.rawGonzDat #storing data into ac
    gonz.calcAK()
