import ipdb
from readData import readUKIRTfields
from readData import readAltbandCrossCheck
import numpy as np
import glob

if __name__=='__main__':

    ipdb.set_trace()
    year = '2017'
    all_files = glob.glob('../data/ukirt/2017/altBandCrossCheck/*'+year+'*')
    for filename in all_files:
        sp_string = filename.split('_')
        field = sp_string[2]+'_'+sp_string[3]
        ccd = sp_string[4].split('.')[0]
        phot = sp_string[1].split('2')[0]
        data = readAltbandCrossCheck(year,field,ccd,photomType=phot,dir='../data/ukirt/2017/altBandCrossCheck/') #delete dir if you need to run Samson
        #make sure value after filename.split('/')[#] (5 right now) is equal to the # of directories
        #in the beginning of newfilename
        newfilename = '../data/ukirt/2017/psfpickles/'+filename.split('/')[5].split('.')[0]+'.npy'
        np.save(newfilename,data)
