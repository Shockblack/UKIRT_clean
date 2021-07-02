from os import path
import numpy as np
#import sys
import plotCMD
#import cmd from plotCMD
import lmfit
import ipdb
import matplotlib.pyplot as plt
from astropy import coordinates as coord
from astropy import units as u
from scipy.signal import find_peaks
import mapper
import parameters as pram
sqrt2pi = np.sqrt(2*np.pi)

class redclumpfinder():
    #def __init__(self,filename,RCnumber,VMIclump,clump,begin=1.5,final=1.5,VMIwidth=0.25,binnumber=100,Iterations=15000):
    #    self.filename=filename
    #def __init__(self,cmdclass,RCnumber,VMIclump,clump,begin=1.5,final=1.5,VMIwidth=0.25,binnumber=100,Iterations=15000):
    def __init__(self,cmdclass,RCnumber=1,VMIclump=1.,clump=999,begin=1.5,final=1.5,VMIwidth=0.25,binnumber=100,Iterations=15000):
        self.cmd = cmdclass
        self.RCnumber=RCnumber
        self.VMIclump=VMIclump
        self.clump=clump
        self.begin=begin
        self.final=final
        self.VMIwidth=VMIwidth
        self.binnumber=binnumber
        self.Iterations=Iterations
        return 

    def fitRCnew(self,func,plotfit=False,A=None,B=None,N_RC=None,M_RC=None,sigma_RC=None):

        #first, find magnitude of RC
        self.fitRCMagnitude(func,plotfit=plotfit,A=A,B=B,N_RC=N_RC,M_RC=M_RC,sigma_RC=sigma_RC)

        A = self.fit.best_values['A']
        B = self.fit.best_values['B']
        M_RC = self.fit.best_values['M_RC']
        sigma_RC = self.fit.best_values['sigma_RC']
        N_RC = self.fit.best_values['N_RC']

        #Next, calculate the weights using this best fit
        self.calcWeights(A,B,M_RC,sigma_RC,N_RC)

        #then, try and determine color
        #bin colors for initial guess
        colorbins = np.linspace(np.percentile(self.cmd.filterStarDict['delta'],5),np.percentile(self.cmd.filterStarDict['delta'],95),50)
        bin_centers = (colorbins[1:]+colorbins[:-1])/2.
        colorhist,colorbins = np.histogram(self.cmd.filterStarDict['delta'],colorbins)
        color_fit_vals = self.determineColor(bin_centers[colorhist.argmax()])
        
        return color_fit_vals


    def determineColor(self,ColorRCinput):
        fitMags = []
        fitHMK = []

        MU1min = max(ColorRCinput-2.5,0)
        MU1max = ColorRCinput-.01
        MU1step = 0.01
        MU1s = np.arange(MU1min,MU1max,step=MU1step)
        MU2 = ColorRCinput
        Kmags = self.cmd.filterStarDict['altmag']
        HMKdiffer = self.cmd.filterStarDict['delta']

        VarianceMin = 1000000

        for MU1 in MU1s:
            
            #calculate differences between of H-K color and foreground disk guess and RCcolorinput
            MU1diff = np.abs(HMKdiffer-MU1)
            MU2diff = np.abs(HMKdiffer-ColorRCinput)
            #find indices that are closer to colorRCinput
            #RCCinds are RCColor indices where MU2diff is less than MU1
            #actually, I'm going to carry these inds into future calcs rather than make new arrays
            #THESE ARE TRUE IF THE COLOR IS CLOSER TO THE RC
            RCCinds = MU1diff>MU2diff
            #those colors closer to RC guess
            HMKverygood = HMKdiffer[RCCinds]
            #weights of those colors closer to RC
            #XXX still need to calc weights
            #weightsverygood = weights[RCCinds]
            #squared color differneces 
            sum2b = HMKdiffer[RCCinds]**2
            #squared sum of color distances for those closer to the foreground disk
            sum2a = HMKdiffer[~RCCinds]**2
            #color status for fitting later, whether rejected or not
            ColorStatus = np.ones(len(HMKdiffer),dtype=bool)
            
            #loop controller
            passer = 1
            passer_counter = 0
            while passer==1 and passer_counter<100:
                passer = 0
                ColorRC = 0
                sumRC = 0
                numberscolor = 0
                combinds = np.logical_and(ColorStatus, RCCinds)
                
                sumRC = np.sum(self.weights[combinds])
                sumColorRC = np.sum(self.weights[combinds]*HMKdiffer[combinds])
                ColorRC = sumColorRC/sumRC
                
                #loop through all the color measurements
                #I think this can be done by ANDing the ColorStatus and RCCinds
                #so only include HMdiffer elements if they pass both

                V1 = np.sum(self.weights[combinds])
                V2 = np.sum(self.weights[combinds]**2)
                VARsum = np.sum(self.weights[combinds]*(HMKdiffer[combinds]-ColorRC)**2)

                #now, get the sigma of the color, unsure, this is variance
                STD2 = (V1/(V1*V1 - V2))*VARsum
                #the cutoff is 2.5sigma
                limitcutoff = 2.5*np.sqrt(STD2)
                #this is sigma
                STDcolor = np.sqrt(STD2)
                #number of failing stars, does nothing
                failingstars = 0
            
                #now we want to copy our combined indices to check if they change later
                preinds = combinds.copy()

            
                #make boolean array for excluding sigma clipped points
                #sigmainds are True for those that pass the sigma clipping
                sigmainds = np.abs(HMKdiffer - ColorRC) < limitcutoff
            
                #ColorStatus is all true going into this
                #this is the only place we modify in the 'passer' loop
                #so, we can check if the sigma clipped point change using colorstatus, and update colorstatus
                #those point which failed before will fail after, and as we narrow through sigma clipping
                if np.array_equal(ColorStatus, np.logical_and(sigmainds, combinds)):
                    #we will AND combine the sigma clipped indices with the colorstatus
                    ColorStatus = np.logical_and(ColorStatus, sigmainds)
                    #also reset passer to be one so we go through the loop again
                    passer = 1
                    passer_counter+=1

                    
                if passer_counter>95:
                    ipdb.set_trace()
            if passer_counter>95:
                #XXX
                #ipdb.set_trace()
                print('passer_counter limit reached')
                return [float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]
                #add logging
            #thus ends passer loop, back in MU1iter
            #printf (" $MU2iter $ColorRC $STDcolor $sumRC\n");
            #factor of 1.096? seems familiar, can't recall why its scaling it up
            varianceB = (1.096*STDcolor)*(1.096*STDcolor);
            RealSTDcolor = 1.096*STDcolor;

            #add up the lower ninty percent of the sum2a array. If its greater than 1, calc variance 
            #otherwise just use the varaince as the calc variance
            sumsquares = np.sum(np.sort(sum2a)[:int(.9*len(sum2a))])
            #foreach(0..int(0.90*$#sum2sort)){
            #for i in 0.9*np.arange(len(sum2sort)):
            #sumsquares += sum2sort[i]
            #numA += 1


            numA = int(.9*len(sum2a))
            if numA > 1:
                varianceA = sumsquares/(numA - 1)*(1.6043)
                VarianceTotal = ((V1 - 1)*varianceB + (numA - 1)*varianceA)/(V1+numA-2)        
            else:
                VarianceTotal  = varianceB
                
            if (VarianceTotal < VarianceMin):
                VarianceMin = VarianceTotal
                MU1opt = MU1
                MU2opt = MU2
                STDtotal = np.sqrt(VarianceTotal)
                FinalColorRC = ColorRC
                RealSTDcolorOpt = RealSTDcolor
                #print('VarianceMin\tSTDTotal\tMU1\tMU2\tRealSTDcolor\tFinalColorRC\n%.3e\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f'%(VarianceMin,STDtotal,MU1opt,MU2opt,RealSTDcolorOpt,FinalColorRC))

        if False: #RealSTDcolorOpt<.4:

            ipdb.set_trace()
            #guesses = np.array(guess_list)
            self.cmd.plotCMD()
            plt.axvline(MU1opt)
            plt.axvline(MU2opt)
            plt.plot(ColorRCinput,self.fit.best_values['M_RC'],'C3o')
            plt.plot(FinalColorRC,self.fit.best_values['M_RC'],'C2o')
            plt.plot([FinalColorRC-RealSTDcolorOpt,FinalColorRC+RealSTDcolorOpt],[self.fit.best_values['M_RC'],self.fit.best_values['M_RC']])
            plt.show()
                                                            
        try:
            return [FinalColorRC, RealSTDcolorOpt, STDtotal, VarianceMin, MU1opt, MU2opt]
        except:
            #XXX VarianceTotal can be nan
            #ipdb.set_trace()
            return [float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]

        #ipdb.set_trace()
        ##Kclump = (
        #Kclump = (M0opt*NRCopt + M0opt2*NRCopt2)/(NRCopt+NRCopt2)
        #if ((np.abs(FinalColorRC - HMKclump) <= 0.03) and (np.abs(Kclump - clump) <= 0.03 )):
        #    HMKclump = FinalColorRC
        #else:
        #    printf ("Restarting!\n");
        
        #printf ("debug here below?\n");
        #printf ("$FinalColorRC \n");
        #printf ("$ColorRC \n");
        #$FinalColorRC = 0.2;
        #printf ("$Iclump \n");
        #system ("perl RedClumpFinder.plx  $filename $RCnumber $FinalColorRC $Iclump $begin $final $VMIwidth $binnumber $Iterations");
        #need to restart
        #exit    
        #ipdb.set_trace()
        #return


    def calcWeights(self,A,B,M_RC,sigma_RC,N_RC):
        #Calculate weights using Nataf 2013 equation 6
        M = self.cmd.filterStarDict['altmag']
        self.weights = RCgaussian(M,M_RC,sigma_RC,N_RC)/redclumpexp(M,A,B,M_RC,sigma_RC,N_RC)
        
        

    def getFitStars(self,starDict):
        #given a cmd class, select stars for later fitting. 
        #filters on color and magnitude cuts, quality cuts are assumed to have happened previously
        #
        #star dict is a dictionary of quality filtered stars
        return

    def icMethod(self,method='mag',numbins=None):
        
        #Calculates the initial conditions for the fitting algorithm
        l = self.cmd.l
        b = self.cmd.b

        #Mag method uses the magnitude and data as IC
        if method=='mag':
            if type(numbins)==type(None):
                numbins = max([20,int(5*np.log(len(self.cmd.fitStarDict['altmag'])))]) #finds number of bins with a minimum of 20
            bins = np.linspace(11,18,numbins)

            val,mbin = np.histogram(self.cmd.fitStarDict['altmag'],bins)
            fitinds = val != 0
            yvals = val[fitinds]
            xvals = ((mbin[1:]+mbin[:-1])/2)[fitinds]

            if l > 350:
                l-=360
            #Using location to aid initial guesses: order here goes outer, inner, midrange (from observation)

            #Outskirts
            if abs(b)>0.9:
                peaks = find_peaks(yvals[np.where(xvals<15)[0]],height=15,width=[8,18])
                if peaks[0].size==0:
                    M_RCguess=13
                else:
                    M_RCguess = xvals[peaks[0][0]]

            #most extinct area close to GC
            elif abs(b) < 0.65 and l < 1.5:
                peaks = find_peaks(yvals,height=70,distance=50)
                if peaks[0].size==0:
                    M_RCguess=15
                else:
                    M_RCguess = xvals[peaks[0][0]]
            #Area in between center and outskirts
            else:
                peaks = find_peaks(yvals,height=70,distance=10)
                if peaks[0].size==0:
                    M_RCguess=14.5
                else:
                    M_RCguess = xvals[peaks[0][0]]
        else:
            exit('Invalid method given')
        return M_RCguess

    def fitRCMagnitude(self,func,plotfit=False,figdir='../misc_figs/maghist2.pdf',A=None,B=None,N_RC=None,M_RC=None,sigma_RC=None,pixID=None):
        
        self.modelMag = lmfit.Model(func)
        fitinds = self.cmd.fitMagHist != 0
        yvals = self.cmd.fitMagHist[fitinds]
        xvals = ((self.cmd.magBins[1:]+self.cmd.magBins[:-1])/2)[fitinds]

        #M_RC = self.icMethod()
        #default guesses from average map parameters
        if type(A)==type(None) or np.isnan(A):# or A<0.2:
            A = 4
        if type(B)==type(None) or np.isnan(B):
            B = .8
        if type(N_RC)==type(None) or np.isnan(N_RC):
            N_RC=yvals.max()
        if type(M_RC)==type(None) or np.isnan(M_RC):
            M_RC = xvals[yvals[np.where(xvals<14.5)[0]].argmax()]            
        if type(sigma_RC)==type(None) or np.isnan(sigma_RC):
            sigma_RC=0.5

        params = self.modelMag.make_params(A=A,B=B,N_RC=N_RC,M_RC = M_RC,sigma_RC=sigma_RC)


        params['A'].min = 0
        params['B'].min = 0
        params['N_RC'].min = 1
        params['M_RC'].min = 12.5
        params['M_RC'].max = 17.
        params['sigma_RC'].min=0
        #params['sigma_RC'].min=.3


        self.fit = self.modelMag.fit(yvals,params,M=xvals)#,weights=1/(yvals))

        if plotfit:
            fig, ax = plt.subplots()

            #test,bins = np.histogram(self.cmd.filterStarDict['altmag']-1*self.cmd.filterStarDict['delta'],self.cmd.magBins)
            #ipdb.set_trace()
            plt.plot(xvals,self.fit.init_fit,'C2--')
            #plt.plot(xvals,yvals,'C0')
            ax.errorbar(xvals,yvals,np.sqrt(yvals),color='C0',linestyle='none')
            #ax.plot(xvals,test[fitinds],color='C2')
            if self.cmd.l > 350:
                self.cmd.l -= 360
            ax.plot(xvals,self.fit.best_fit,'C1--',linewidth=3)
            if type(pixID) != type(None): #For if I want to know know what pixel this specific plot comes from. Useful for debugging
                ax.set_title('(l,b)=(%.4f,%.4f), size=%0.1f$^\prime$, Pix=%d'%(self.cmd.l,self.cmd.b,60*self.cmd.edge_length,pixID))
            else:
                ax.set_title('(l,b)=(%.4f,%.4f), size=%0.1f$^\prime$'%(self.cmd.l,self.cmd.b,60*self.cmd.edge_length))
            ax.axvline(self.fit.best_values['M_RC'],color='k',linestyle='--',linewidth=3)
            ax.set_ylabel('Number')
            ax.set_xlabel('$\it{K}$ Magnitude')
            ax.text(self.fit.best_values['M_RC']+.3,5,'$\it{K_{RC}}=$%.2f'%(self.fit.best_values['M_RC']),fontsize=16)
            plt.savefig(figdir)
            #plt.show()

        return self.fit.best_fit


    def fitRCcolor(self):
        return


def RCgaussian(M,M_RC,sigma_RC,N_RC):
    #calculate just the RC gaussian from fitted parameters
    return N_RC/(sqrt2pi*sigma_RC)*np.exp(-(M-M_RC)**2/(2*sigma_RC**2))


def redclumpexp(M,A,B,M_RC,sigma_RC,N_RC):
        
        #definitions

    N_RGBB = 0.201*N_RC
    N_AGBB = 0.028*N_RC
    M_RGBB = M_RC+0.737
    M_AGBB = M_RC-1.07
    sigma_RGBB = sigma_RC
    sigma_AGBB = sigma_RC
    
    term1 = A*np.exp(B*(M-M_RC))
        #Gaussian for the red clump
    term2 = N_RC/(sqrt2pi*sigma_RC)*np.exp(-(M-M_RC)**2/(2*sigma_RC**2))
        #Gaussian for red giant branch bump
    term3 = N_RGBB/(sqrt2pi*sigma_RGBB)*np.exp(-(M-M_RGBB)**2/(2*sigma_RGBB**2))
        #Gaussian for asymptotic giant branch bump
    #term4 = N_AGBB/(sqrt2pi*sigma_AGBB)*np.exp(-(M-M_AGBB)**2/(2*sigma_AGBB**2))
    
    NMdM = term1 + term2 + term3# + term4
    return NMdM

def redclumpOnlyExp(M,A,B,M_RC,sigma_RC,N_RC):
        
        #definitions

    N_RGBB = 0.201*N_RC
    N_AGBB = 0.028*N_RC
    M_RGBB = M_RC+0.737
    M_AGBB = M_RC-1.07
    sigma_RGBB = sigma_RC
    sigma_AGBB = sigma_RC
    
    term1 = A*np.exp(B*(M-M_RC))
        #Gaussian for the red clump
    term2 = N_RC/(sqrt2pi*sigma_RC)*np.exp(-(M-M_RC)**2/(2*sigma_RC**2))
        #Gaussian for red giant branch bump
    term3 = N_RGBB/(sqrt2pi*sigma_RGBB)*np.exp(-(M-M_RGBB)**2/(2*sigma_RGBB**2))
        #Gaussian for asymptotic giant branch bump
    #term4 = N_AGBB/(sqrt2pi*sigma_AGBB)*np.exp(-(M-M_AGBB)**2/(2*sigma_AGBB**2))
    
    NMdM = term1 + term2 + term3# + term4
    return NMdM
    
#def redclumpclassic(self,M,a,b,c,M_RC,sigma_RC):
def redclumpclassic(M,a,b,c,N_RC,M_RC,sigma_RC):
    quad = a + b*(M-M_RC) + c*(M-M_RC)**2
    gauss = N_RC/(sqrt2pi*sigma_RC) * np.exp(-(M-M_RC)**2/(2*sigma_RC**2))
    NMdM = quad+gauss
    return NMdM



if __name__=='__main__':
    
    #ipdb.set_trace()
    #cmd_test = plotCMD.cmd('test.txt',268.5,-28.7,edge_length=3/60.)
    #cmd_test = plotCMD.cmd('test.txt',269.385168,-28.971528,edge_length=2/60.)
    #cmd_test = plotCMD.cmd('test.txt',269.318501,-29.771528,edge_length=2/60.)
    #cmd_test = plotCMD.cmd('test.txt',268.985168,-29.604862,edge_length=2/60.)
    #cmd_test = plotCMD.cmd('test.txt',268.918501,-29.738195,edge_length=2/60.)
    #cmd_test = plotCMD.cmd('test.txt',268.885168,-29.704862,edge_length=2/60.)
    #cmd_test = plotCMD.cmd('test.txt',269.36695 , -28.98901,edge_length=3/60.)
    #cmd_test = plotCMD.cmd('test.txt',269.36766 , -29.05988,edge_length=3/60.)
    #cmd_test = plotCMD.cmd('test.txt',266.418501 , -29.004862,edge_length=2/60.)
    #cmd_test = plotCMD.cmd('test.txt',t,edge_length=2/60.)
    #cmd_test = plotCMD.cmd('test.txt',266.,-29.,edge_length=3/60.)
    #cmd_test = plotCMD.cmd('test.txt',266.968501,-29.038195,edge_length=2/60.)
    cmd_test = plotCMD.cmd('test.txt',266.,-29.,edge_length=pram.arcmin/60.)
    #cmd_test = plotCMD.cmd('test.txt',269.,-29.3,edge_length=.1)
    #cmd_test = plotCMD.cmd('test.txt',268.5,-29.7,edge_length=.1)

    rcfinder = redclumpfinder(cmd_test)
    M_RCguess = rcfinder.icMethod()
    rcfinder.fitRCnew(redclumpOnlyExp,plotfit=False,M_RC=M_RCguess)
    #ipdb.set_trace()
    rcfinder.fitRCMagnitude(redclumpOnlyExp, plotfit=True, M_RC=M_RCguess)

    #ipdb.set_trace()

    rcfinder.cmd.color_fit_cut(rcfinder.fit.best_values['M_RC'],rcfinder.fit.best_values['sigma_RC'],2)
    #rcfinder.cmd.plotCMD()
    print(rcfinder.fit.fit_report())
    #ipdb.set_trace()
    print('Number of RC stars in fit: %f\n'%(np.sqrt(np.pi/(2*rcfinder.fit.best_values['sigma_RC']**2))*rcfinder.fit.best_values['N_RC']/(sqrt2pi*rcfinder.fit.best_values['sigma_RC'])))
    #ipdb.set_trace()
    """
    if len(sys.argv) == 4:
        filename = $ARGV[0]
        RCnumber = $ARGV[1]
        VMIclump = $ARGV[2]
        clump = $ARGV[3]
        begin = 1.5
        final = 1.5
        VMIwidth = 0.25
        binnumber = 100
        Iterations = 15000

    elif len(sys.argv) == 8 :
        filename = $ARGV[0]
        RCnumber = $ARGV[1]
        VMIclump = $ARGV[2]
        clump = $ARGV[3]
        begin = $ARGV[4]
        final = $ARGV[5]
        VMIwidth = $ARGV[6]
        binnumber = $ARGV[7]
        Iterations = $ARGV[8]
    """
