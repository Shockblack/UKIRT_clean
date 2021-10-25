from operator import le
import numpy as np
from astropy import coordinates as coord
from astropy import units as u
from readData import readUKIRTfields
import ipdb
import RedClumpFinder as rcf
import plotCMD as pcmd
import random as ran
import parameters as pram

class map:
    def __init__(self,edge_length=pram.arcmin/60.,year=pram.year,ra_lims=None,dec_lims=None):
        #limits in case we want to make a smaller map
        self.edge_length = edge_length
        self.ra_lims = ra_lims
        self.dec_lims = dec_lims
        self.year = year
        #ipdb.set_trace()
        self.raw_fieldData = readUKIRTfields()
        
        #select by year, probably shouldn't go here
        self.fieldData = []
        for field in self.raw_fieldData:
            if int(field['year']) == year:
                self.fieldData.append(field)

    #function to read total field data
    def get_fields_stats(self):
        if (self.ra_lims != None) and (self.dec_lims != None):
            self.ra_min = self.ra_lims[0]
            self.ra_max = self.ra_lims[1]

            self.dec_min = self.dec_lims[0]
            self.dec_max = self.dec_lims[1]
            return
            print('test')
        
        self.allRAmax = []
        self.allRAmin = []
        self.allDECmax = []
        self.allDECmin = []
        #ipdb.set_trace()
        for field in self.fieldData:
            self.allRAmax.append(field['RAmax'])
            self.allRAmin.append(field['RAmin'])
            self.allDECmax.append(field['DECmax'])
            self.allDECmin.append(field['DECmin'])
        
        self.allRAmax = np.array(self.allRAmax)
        self.allRAmin = np.array(self.allRAmin)
        self.allDECmax = np.array(self.allDECmax)
        self.allDECmin = np.array(self.allDECmin)

        self.ra_min = self.allRAmin.min()
        self.ra_max = self.allRAmax.max()

        self.dec_min = self.allDECmin.min()
        self.dec_max = self.allDECmax.max()

        return

    #function to generate pixel grid
    def gen_grid(self,grid_grow_scale=0.025):
        
        #ipdb.set_trace()
        ra_span = self.ra_max-self.ra_min
        dec_span = self.dec_max-self.dec_min

        ra_max = self.ra_max+grid_grow_scale*ra_span
        ra_min = self.ra_min-grid_grow_scale*ra_span

        dec_max = self.dec_max+grid_grow_scale*dec_span
        dec_min = self.dec_min-grid_grow_scale*dec_span

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

    
    def filter_grid(self):
        npixels = 0
        self.pixel_centers = []
        for pixel in self.grid_pixel_centers:
            for field in self.fieldData:
                #if the pixel falls in a field, add it to the list
                if pixel[0]<=field['RAmax'] and pixel[0]>=field['RAmin'] and pixel[1]<=field['DECmax'] and pixel[1]>=field['DECmin']:
                    self.pixel_centers.append(pixel)
                    npixels+=1
                    #pixel passsed, no more checking needed
                    break
        self.npixels = npixels
        print('Got '+str(npixels)+' to fit!')


    def fit_map(self,filebase='test_map',method='pix',plotmag=False,checkpix=True):
        #ipdb.set_trace()
        self.pixels = []

        #Following lines are variables used in calculating mag difference in
        #lines on map to show UKIRT team
        lowerMagAvg = 0
        lowerSum = 0
        upperMagAvg = 0
        upperSum = 0
        normMagAvg = 0
        normSum = 0
        ilist = []

        #first pixel done manually
        pixel = self.pixel_centers[0]
        cmd = pcmd.cmd('test.txt',pixel[0],pixel[1],l=pixel[2],b=pixel[3],edge_length=self.edge_length)

        #ratio used as a counter for plotting mag histograms
        ratio = ran.randint(1,int(self.npixels/30))
        badpix = []
        i=0
        #for all but first pixel
        for pixel in self.pixel_centers[::-1]:
            i+=1
            self.pixels.append(pixel)
            self.pixels[-1].append(self.edge_length)
            #ipdb.set_trace()
            if cmd.check_new_coords(pixel[0],pixel[1]):
                cmd.updateCMD(pixel[0],pixel[1],pixel[2],pixel[3])
            else:
                del cmd
                cmd = pcmd.cmd('test.txt',pixel[0],pixel[1],l=pixel[2],b=pixel[3],edge_length=self.edge_length)
            #if i == 10:
            #    ipdb.set_trace()
            rcfinder = rcf.redclumpfinder(cmd)
            self.pixels[-1].append(len(rcfinder.cmd.fitStarDict['name']))
            #rcfinder.fitRCMagnitude(rcf.redclumpOnlyExp,plotfit=False)
            #header = 'RA,DEC,edgelength,N_stars,A,Aerr,B,Berr,M_RC,M_RCerr,sigma_RC,sigma_RCerr,N_RC,N_RCerr,FinalColorRC,RealSTDcolorOpt,STDtotal,VarianceMin,MU1opt,MU2opt'
            print('Fitting')
            #============================================================================
            #Pix method is an old method that determines the initial conditions based on
            #those of the previous pixel. Good idea in theory, but didn't work out well
            #overall. Currently use the actual data from the CMD to determine initial conditions.
            #============================================================================
            #our two options
            #method='pix'

            if pixel[2]>358.96 and pixel[2] < 359.0113 and pixel[3] > -1.52 and pixel[3] < -1.436:
                ilist.append(i)

            method='mag'
            if method=='mag':

                #numbins = max([30,int(5*np.log(len(rcfinder.cmd.fitStarDict['altmag'])))])
                numbins = 90
                relunc = 1
                relunc1 = 0.5
                numbins1 = numbins

                M_RCguess = rcfinder.icMethod(numbins=numbins)
                color_fit_vals = rcfinder.fitRCnew(rcf.redclumpOnlyExp,M_RC=M_RCguess)

#                while relunc > 0.25 and numbins > 20: #Changes amount of bins depending on the uncertainty 
#                    #IC are now calculated in RedClumpFinder.py, which the line below does
#                    M_RCguess = rcfinder.icMethod(numbins=numbins)
#                    
#                    color_fit_vals = rcfinder.fitRCnew(rcf.redclumpOnlyExp,M_RC=M_RCguess)
#
#                    try:
#                        relunc = np.sqrt(rcfinder.fit.covar[2,2])/rcfinder.fit.best_values['M_RC']
#                    except:
#                        pass
#                    numbins -= 10
#
#                    #This section is just making sure that the unc is indeed decreasing with bin decrease
#                    M_RCguess = rcfinder.icMethod(numbins=numbins)
#                    color_fit_vals = rcfinder.fitRCnew(rcf.redclumpOnlyExp,M_RC=M_RCguess)
#                    try:
#                        relunc1 = np.sqrt(rcfinder.fit.covar[2,2])/rcfinder.fit.best_values['M_RC']
#                        if relunc1 > relunc or relunc1 == relunc:
#                            print("Uncertainty did not decrease!")
#                    except:
#                        pass
#                    
                while relunc1 < relunc and numbins >20:
                    numbins = numbins1
                    M_RCguess = rcfinder.icMethod(numbins=numbins)
                    color_fit_vals = rcfinder.fitRCnew(rcf.redclumpOnlyExp,M_RC=M_RCguess)
                    try:
                        relunc = np.sqrt(rcfinder.fit.covar[2,2])/rcfinder.fit.best_values['M_RC']
                    except:
                        pass
                    
                    numbins1 = numbins-5

                    M_RCguess = rcfinder.icMethod(numbins=numbins1)
                    color_fit_vals = rcfinder.fitRCnew(rcf.redclumpOnlyExp,M_RC=M_RCguess)
                    try:
                        relunc1 = np.sqrt(rcfinder.fit.covar[2,2])/rcfinder.fit.best_values['M_RC']
                        if relunc1 > relunc or relunc1 == relunc:
                            print("Uncertainty did not decrease!")
                    except:
                        pass

#                try:
#                    rchisq = rcfinder.fit.redchi
#                except:
#                    print('NO RCHISQR')
#
#                numbins1=numbins-5
#                M_RCguess = rcfinder.icMethod(numbins=numbins1)
#                color_fit_vals = rcfinder.fitRCnew(rcf.redclumpOnlyExp,M_RC=M_RCguess)
#                try:
#                    rchisq1 = rcfinder.fit.redchi
#                except:
#                    print('NO RCHISQR1')
#
#                while abs(rchisq1-1)<abs(rchisq-1) and numbins > 20:
#                    numbins=numbins1
#                    M_RCguess = rcfinder.icMethod(numbins=numbins1)
#                    color_fit_vals = rcfinder.fitRCnew(rcf.redclumpOnlyExp,M_RC=M_RCguess)
#                    try:
#                        rchisq = rcfinder.fit.redchi
#                    except:
#                        pass
#
#                    numbins1 = numbins - 5
#                    M_RCguess = rcfinder.icMethod(numbins=numbins1)
#                    color_fit_vals = rcfinder.fitRCnew(rcf.redclumpOnlyExp,M_RC=M_RCguess)
#
#                    try:
#                        rchisq1 = rcfinder.fit.redchi 
#                        if abs(rchisq1-1)>abs(rchisq-1) or rchisq == rchisq1:
#                            print("REDCHI did not decrease!")
#                    except:
#                        pass
                
                if plotmag == True:
                    if i == ratio:
                        rcfinder.fitRCMagnitude(rcf.redclumpOnlyExp,plotfit=True,figdir='../misc_figs/maghist'+'_'+str(i)+'.pdf',M_RC=M_RCguess)
                        ratio+=int(self.npixels/30)
                
#                if relunc > 0.25:
#                    rcfinder.fitRCMagnitude(rcf.redclumpOnlyExp,plotfit=True,figdir='../misc_figs/maghist'+'_'+str(i)+'.pdf',M_RC=M_RCguess)
                print('Done fitting!')
                
            elif method=='pix':
                try:
                    
                    #check the geometric distance of the last pixel to make sure it is right next to the new pixel
                    if np.sqrt((self.pixels[-1][0]-self.pixels[-2][0])**2+(self.pixels[-1][1]-self.pixels[-2][1])**2) < 3*self.edge_length and self.pixels[-2][13]/self.pixels[-2][12]<.008:
                        print('Using last pixel as init')
                        color_fit_vals = rcfinder.fitRCnew(rcf.redclumpOnlyExp,A=self.pixels[-2][6],B=self.pixels[-2][8],N_RC=self.pixels[-2][14],M_RC=self.pixels[-2][10],sigma_RC=self.pixels[-2][12])
                    else:
                        #Defaults to average values in RCF
                        color_fit_vals = rcfinder.fitRCnew(rcf.redclumpOnlyExp)            
                except:
                    #Defaults to average values in RCF
                    print('Init guess failed')
                    color_fit_vals = rcfinder.fitRCnew(rcf.redclumpOnlyExp)            
                print('Done fitting!')
                
            else:
                exit('Invalid method given')
            
            
            #self.pixels[-1].append(np.sqrt(self.pixels[-1][0]**2+self.pixels[-1][1]**2))
            j=0
            for key in rcfinder.fit.best_values.keys():
                self.pixels[-1].append(rcfinder.fit.best_values[key])
                try:
                    self.pixels[-1].append(np.sqrt(rcfinder.fit.covar[j,j]))
                except:
                    self.pixels[-1].append(float('NaN'))
                j+=1
                #ipdb.set_trace()
            for val in color_fit_vals:
                self.pixels[-1].append(val)
            try:
                #if checkpix == True and abs(self.pixels[-1][10]-self.pixels[-2][10]) > 1.2 and abs(self.pixels[-1][10]-self.pixels[-3][10]) > 1.2 and \
                if checkpix == True and self.pixels[-1][10] < 14 and abs(self.pixels[-1][3]) < 0.3 and self.pixels[-1][2] <0.5 and \
                    np.sqrt((self.pixels[-1][0]-self.pixels[-2][0])**2+(self.pixels[-1][1]-self.pixels[-2][1])**2) < 3*self.edge_length:
                    pixcheck = []
                    rcfinder.fitRCMagnitude(rcf.redclumpOnlyExp,plotfit=True,figdir='../misc_figs/BADmaghist'+'_'+str(i)+'.pdf',M_RC=M_RCguess,pixID=i)
                    pixcheck.append(i)
                    pixcheck.append(pixel[0])
                    pixcheck.append(pixel[1])
                    pixcheck.append(pixel[2])
                    pixcheck.append(pixel[3])
                    badpix.append(pixcheck)
            except:
                pass
            print('pixel '+str(i)+' of '+str(self.npixels)+' done!')

            checkarea = 3
            l = pixel[2]
            if l > 350:
                l -= 360

            if i >= 4637 and i <= 4665:
                normSum += self.pixels[-1][10]
                normMagAvg = normSum/(i-4739)

            if i >= 4740 and i <= 4768:
                lowerSum += self.pixels[-1][10]
                lowerMagAvg = lowerSum/(i-4739)

            if i >= 4848 and i <= 4876:
                upperSum += self.pixels[-1][10]
                upperMagAvg = upperSum/(i-4847)

            if checkarea == 1:
            #Creates a cmd with histogram subplot if the number of red clump stars is too low
                if self.pixels[-1][14] < 10 and abs(pixel[3]) < 0.65 and l < 1.5:
                    limit_dict = {'N':[10,1e6],'altN':[3,1000]}# {'mag':[12,15]}
                    cmd.getStars(limit_dict)
                    cm_dict = {'altmag':[12,16],'delta':[-1,5]}
                    cmd.color_mag_cut(cm_dict,percentile_cut=True)
                    rcfinder=rcf.redclumpfinder(cmd)
                    M_RCguess = rcfinder.icMethod()
                    fit = rcfinder.fitRCMagnitude(rcf.redclumpOnlyExp, plotfit=False, M_RC=M_RCguess)
                    cmd.plotCMDhist(cm_dict,plotsave=True,fit=fit,figdir='../misc_figs/BadNRC/badNRC_'+str(i)+'.pdf',pixID=i)
            elif checkarea == 2:
                if self.pixels[-1][10] > 14.7 and pixel[3] < -0.8 and pixel[3] > -1.4 and l > 0.8:
                    limit_dict = {'N':[10,1e6],'altN':[3,1000]}# {'mag':[12,15]}
                    cmd.getStars(limit_dict)
                    cm_dict = {'altmag':[12,16],'delta':[-1,5]}
                    cmd.color_mag_cut(cm_dict,percentile_cut=True)
                    rcfinder=rcf.redclumpfinder(cmd)
                    M_RCguess = rcfinder.icMethod()
                    fit = rcfinder.fitRCMagnitude(rcf.redclumpOnlyExp, plotfit=False, M_RC=M_RCguess)
                    cmd.plotCMDhist(cm_dict,plotsave=True,fit=fit,figdir='../misc_figs/BadNRC/badNRC_'+str(i)+'.pdf',pixID=i)

            elif checkarea == 3:
                #if self.pixels[-1][10] > 14.7 and l > 1.5:
                if i ==6730:
                    limit_dict = {'N':[10,1e6],'altN':[3,1000]}# {'mag':[12,15]}
                    cmd.getStars(limit_dict)
                    cm_dict = {'altmag':[12,16],'delta':[-1,5]}
                    cmd.color_mag_cut(cm_dict,percentile_cut=True)
                    rcfinder=rcf.redclumpfinder(cmd)
                    M_RCguess = rcfinder.icMethod()
                    fit = rcfinder.fitRCMagnitude(rcf.redclumpOnlyExp, plotfit=False, M_RC=M_RCguess)
                    cmd.plotCMDhist(cm_dict,plotsave=True,fit=fit,figdir='../misc_figs/BadNRC/badNRC_'+str(i)+'.jpg')#,pixID=i)
            #del rcfinder
            #ipdb.set_trace()
            if i%1000 == 0:
                self.save_map(filename=filebase+'.'+str(i)+'.tmp')
            
        print('Pixels for review including [Pix ID, RA, DEC, l, b] are:')
        for k in range(len(badpix)):
            print(badpix[k])
        print("Upper mag avg: ",upperMagAvg)
        print("Lower mag avg: ",lowerMagAvg)
        print("Norm mag avg: ", normMagAvg)
        print("Weird Mag Difference: ", upperMagAvg-lowerMagAvg)
        print('Norm Mag Difference: ', abs(lowerMagAvg-normMagAvg))

        print(ilist)
        ipdb.set_trace()



    def fill_map(self):
        #ipdb.set_trace()
        self.pixels = []
        for pixel in self.pixel_centers:
            self.pixels.append(pixel)
            self.pixels[-1].append(self.edge_length)
            self.pixels[-1].append(np.sqrt(self.pixels[-1][0]**2+self.pixels[-1][1]**2))
            for i in range(9):
                self.pixels[-1].append(np.random.uniform())
        #ipdb.set_trace()


    #function to optimize grid position to minimize empty space
    def optimize_grid(self):
        #remove those with overlap != 1? 
        return
    

    #function to save map
    #output will be lines with field center RA, Dec, edgelength, {params}
    def save_map(self,filename='test_map',path='maps/'):
        header = 'RA,DEC,l,b,edgelength,N_stars,A,Aerr,B,Berr,M_RC,M_RCerr,sigma_RC,sigma_RCerr,N_RC,N_RCerr,FinalColorRC,RealSTDcolorOpt,STDtotal,VarianceMin,MU1opt,MU2opt'
        #date 
        #year edge_length number_pixels 
        #header for columns
        #we will probably want some header info
        with open(path+filename+'.map','w') as outfile:
            for pixel in self.pixels:
                write_string = ""
                #want to pre-calc l,b, have those as 2 and 3
                for item in pixel:
                    write_string+='%0.6f,'%item
                #ipdb.set_trace()
                outfile.write(write_string[:-1]+'\n')
        return

if __name__=='__main__':
    #ipdb.set_trace()
    test_map = map(year=pram.year,edge_length=pram.arcmin/60.)#,ra_lims=[269,280],dec_lims=[-31,-26])#,ra_lims=[268.75,289.25],dec_lims=[-29.75,-29.25])
    test_map.get_fields_stats()
    test_map.gen_grid()
    test_map.filter_grid()
    #Converting inputs to strings so file name is automatically adjusted for inputs
    yearstr = str(pram.year)
    aminedge = str(pram.arcmin)
    filename = 'map_'+pram.phot+'_'+yearstr+'_'+aminedge
    #filename = 'inve_guess_2017_2'
    #ipdb.set_trace()
    test_map.fit_map(filebase=filename,plotmag=False,checkpix=True)
    #ipdb.set_trace()
    test_map.save_map(filename=filename)
    #ipdb.set_trace()
