#-----------------------------------------------------------------------
#
# File name: RC_fitter.py
#
# Replacement file for RedClumpFinder.py which hopes to remove unneeded
# functions and allow for a lot more flexibility in the fitting process.
# Main goal is to be able to refit the RC with different parameters and
# after each each.
#
# Programmer: Aiden Zelakiewicz (zelakiewicz.1@osu.edu)
#
# Revision History:
#   12-Oct-2022 :   File Created
#-----------------------------------------------------------------------

# Importing all the necessary packages

from astropy import coordinates as coord
from astropy import units as u
from scipy import integrate, LowLevelCallable, optimize
import matplotlib.pyplot as plt
import numpy as np
import ipdb
import lmfit

import emcee
import corner

import time
import os
import ctypes

import parameters as pram
import createCMD

sqrt2pi = np.sqrt(2*np.pi)

class RedClump():
    def __init__(self, cmd, binnumber=100, nwalkers=200, iterations=1000, burnin=100):
        self.cmd = cmd
        self.binnumber = int(.05*len(cmd.fitStarDict['altmag']))
        self.nwalkers = nwalkers
        self.iterations = iterations
        self.burnin = burnin

        self.N_obs = len(cmd.fitStarDict['altmag'])

        self.M, self.N, self.Nerr = self.prepare_data_hist(self.cmd.fitStarDict['altmag'])

        # Importing the compiled C function to integrate
        self.lib = ctypes.CDLL(os.path.abspath('integratelib.so'))
        self.lib.f.restype = ctypes.c_double
        self.lib.f.argtypes = (ctypes.c_int, ctypes.POINTER(ctypes.c_double), ctypes.c_void_p)


    def model_MCMC(self, theta, M):
        """The model for the MCMC fitting.

        Parameters
        ----------
        theta : array
            Array of the parameters for the model.
        M : array
            Array of the magnitudes of the stars.

        Returns
        -------
        array
            Array of the model for the MCMC fitting.
        """
        A, B, M_RC, sigma_RC, N_RC = theta

        # EW, B, M_RC, sigma_RC = theta

        # A = From integral

        # N_RC = EW*A

        # Move log space to linear

        return luminosityFunction(M,A,B,M_RC,sigma_RC,N_RC)

    def model_MCMC_nataf(self, theta, M):
        """The model for the MCMC fitting.

        Parameters
        ----------
        theta : array
            Array of the parameters for the model.
        M : array
            Array of the magnitudes of the stars.

        Returns
        -------
        array
            Array of the model for the MCMC fitting.
        """
        EWRC, B, M_RC, sigma_RC = theta['EWRC'].value, theta['B'].value, theta['MRC'].value, theta['SIGMA'].value

        A = self.A

        N_RC = EWRC*A

        # Move log space to linear

        return luminosityFunction(M,A,B,M_RC,sigma_RC,N_RC)


    def log_prior(self, theta):
        """The prior for the MCMC fitting.

        Parameters
        ----------
        theta : array
            Array of the parameters for the model.

        Returns
        -------
        float
            The log of the prior.
        """
        A, B, M_RC, sigma_RC, N_RC = theta
        if (1< A < 20 and 0 < B < 2 and 12 < M_RC < 17 and 0 < sigma_RC < 1 and 0 < N_RC < 100):
            return 0.0
        return -np.inf

    def log_prior_nataf(self, theta):
        """The prior for the MCMC fitting.

        Parameters
        ----------
        theta : array
            Array of the parameters for the model.

        Returns
        -------
        float
            The log of the prior.
        """
        EWRC, B, M_RC, sigma_RC = theta
        if (0< EWRC < 4 and 0 < B < 2 and 12 < M_RC < 17 and 0 < sigma_RC < 1):
            return 0.0
        return -np.inf

    def log_likelihood(self, theta, M, N, Nerr):
        """The likelihood for the MCMC fitting.

        Parameters
        ----------
        theta : array
            Array of the parameters for the model.
        M : array
            Array of the magnitudes of the stars.
        N : array
            Array of the number of stars in each magnitude bin.
        Nerr : array
            Array of the error in the number of stars in each magnitude bin.

        Returns
        -------
        float
            The log of the likelihood.
        """
        # import ipdb; ipdb.set_trace()
        model = self.model_MCMC(theta, M)
        Nobs = len(self.cmd.fitStarDict['altmag'])
        # return np.sum(np.log(model))-Nobs
        return -0.5*np.sum(((N-model)/Nerr)**2)

    def log_likelihood_nataf(self, theta, M):
        # import ipdb; ipdb.set_trace()

        EWRC, B, M_RC, sigma_RC = theta['EWRC'].value, theta['B'].value, theta['MRC'].value, theta['SIGMA'].value

        # Preparing the scipy integrator in C
        py_vals = [EWRC, B, M_RC, sigma_RC]

        self.value_tracker.append(py_vals)
        # ipdb.set_trace()
        c_vals = (ctypes.c_double * len(py_vals))(*py_vals)
        user_data = ctypes.cast(ctypes.pointer(c_vals), ctypes.c_void_p)
        func = LowLevelCallable(self.lib.f, user_data)

        I = integrate.quad(func, self.cmd.cm_dict['altmag'][0], self.cmd.cm_dict['altmag'][1])
        if I[0]==0:
            ipdb.set_trace()

        self.A = self.N_obs/I[0]#self.integrator(EWRC, B, M_RC, sigma_RC, n=100)

        model = self.model_MCMC_nataf(theta, M)

        return np.sum(np.log(model))-self.N_obs

    def log_likelihood_nataf_neg(self, theta, M):
        # import ipdb; ipdb.set_trace()

        EWRC, B, M_RC, sigma_RC = theta['EWRC'].value, theta['B'].value, theta['MRC'].value, theta['SIGMA'].value

        # Preparing the scipy integrator in C
        py_vals = [EWRC, B, M_RC, sigma_RC]

        self.value_tracker.append(py_vals)
        # ipdb.set_trace()
        c_vals = (ctypes.c_double * len(py_vals))(*py_vals)
        user_data = ctypes.cast(ctypes.pointer(c_vals), ctypes.c_void_p)
        func = LowLevelCallable(self.lib.f, user_data)

        I = integrate.quad(func, self.cmd.cm_dict['altmag'][0], self.cmd.cm_dict['altmag'][1])
        if I[0]==0:
            ipdb.set_trace()

        self.A = self.N_obs/I[0]#self.integrator(EWRC, B, M_RC, sigma_RC, n=100)

        model = self.model_MCMC_nataf(theta, M)

        return -(np.sum(np.log(model))-self.N_obs)

    def log_probability(self, theta, M):#, N, Nerr):
        """The probability for the MCMC fitting.

        Parameters
        ----------
        theta : array
            Array of the parameters for the model.
        M : array
            Array of the magnitudes of the stars.
        N : array
            Array of the number of stars in each magnitude bin.
        Nerr : array
            Array of the error in the number of stars in each magnitude bin.

        Returns
        -------
        float
            The log of the probability.
        """
        lp = self.log_prior_nataf(theta)
        if not np.isfinite(lp):
            return -np.inf
        return lp + self.log_likelihood_nataf(theta, M)#, N, Nerr)

    def prepare_data_hist(self, M):
        """Prepares the data for the histogram.

        Parameters
        ----------
        M : array
            Array of the magnitudes of stars.

        Returns
        -------
        array
            Array of the magnitudes of the stars.
        array
            Array of the number of stars in each magnitude bin.
        array
            Array of the error estimation in the number of stars in each magnitude bin.
            *THIS IS NOT FINAL ERRORS, JUST ESTIMATES FOR THE MCMC FITTING*
        """
        bins = np.linspace(12,16,self.binnumber)
        Nstars, hist_bins = np.histogram(M,bins)#'fd')

        mags = (hist_bins[1:]+hist_bins[:-1])/2
        Nstars_err_guess = np.sqrt(Nstars)

        Nstars_err_guess[Nstars_err_guess == 0] = 5
        
        return mags, Nstars, Nstars_err_guess

    def run_MCMC(self, M, initial_guess=None):#, N, Nerr):
        """Runs the MCMC fitting.

        Parameters
        ----------
        M : array
            Array of the magnitudes of stars.
        N : array
            Array of the number of stars in each magnitude bin.
        Nerr : array
            Array of the error estimation in the number of stars in each magnitude bin.
            *THIS IS NOT FINAL ERRORS, JUST ESTIMATES FOR THE MCMC FITTING*

        Returns
        -------
        array
            Array of the parameters for the model.
        array
            Array of the covariance matrix for the parameters.
        """
        # Check if we are given an initial guess from Nelder-Mead
        if type(initial_guess) == type(None):
            initial_guess = np.array([1.2,0.8,14,0.5])

        p0 = [np.array(initial_guess) + np.random.randn(len(initial_guess)) for i in range(self.nwalkers)]

        sampler = emcee.EnsembleSampler(self.nwalkers, len(initial_guess), self.log_probability, args=(M, ))

        # Run the burn-in
        p0, _, _ = sampler.run_mcmc(p0, self.burnin)
        sampler.reset()

        # Run the production
        pos, prob, state =  sampler.run_mcmc(p0, self.iterations, progress=True)

        # Get the best fit parameters
        best_ind = np.argmax(sampler.flatlnprobability)
        samples = sampler.flatchain
        best_params = samples[best_ind]

        # Calculate the statistical uncertainties
        perc = np.percentile(samples, [16, 50, 84], axis=0)
        unc = np.diff(perc, axis=0).T
        unc = np.sqrt(unc[:,0]**2 + unc[:,1]**2)

        # Calculating the best integral
        # Preparing the scipy integrator in C
        c_vals = (ctypes.c_double * len(best_params))(*best_params)
        user_data = ctypes.cast(ctypes.pointer(c_vals), ctypes.c_void_p)
        func = LowLevelCallable(self.lib.f, user_data)

        I = integrate.quad(func, self.cmd.cm_dict['altmag'][0], self.cmd.cm_dict['altmag'][1])

        A = self.N_obs/I[0]
        # Propogating the error
        A_err = A*np.sqrt( (np.sqrt(self.N_obs)/self.N_obs)**2 + (I[1]/I[0])**2 )
        N_RC = A*best_params[0]
        N_RC_err = N_RC*np.sqrt((A_err/A)**2 + (unc[0]/best_params[0])**2)

        params =    {'A':A, 'A_err':A_err, 'B':best_params[1], 'B_err':unc[1], 'MRC':best_params[2], 'MRC_err':unc[2],\
                    'SIGMA':best_params[3], 'SIGMA_err':unc[3], 'NRC':N_RC, 'NRC_err':N_RC_err}

        return sampler, params, unc

    def fit(self):
        M, N, Nerr = self.prepare_data_hist(self.cmd.fitStarDict['altmag'])
        sampler, pos, prob, state = self.run_MCMC(M, N, Nerr)
        best_ind = np.argmax(sampler.flatlnprobability)
        samples = sampler.flatchain
        best_params = samples[best_ind]
        return best_params

    def find_color(self):
        return

    def integrator(self, EWRC, B, M_RC, sigma_RC, n=1000):
        """
        Integrates the luminosity function to find the normalization constant A
        using the trapezoidal rule.
        """
        a = self.cmd.cm_dict['altmag'][0]
        b = self.cmd.cm_dict['altmag'][1]
        
        h = (b-a)/n

        I = 0.5*(luminosityFunction_EWRC(a, EWRC, B, M_RC, sigma_RC) + luminosityFunction_EWRC(b, EWRC, B, M_RC, sigma_RC))
        for k in range(1, n):
            I += luminosityFunction_EWRC(a+k*h, EWRC, B, M_RC, sigma_RC)
        I = I*h
        return I

    def run_minimizer(self, M, initial_guess=None):

        self.value_tracker = []

        if type(initial_guess) == type(None):
            initial_guess = np.array([1.2,0.8,14,0.5])

        params = lmfit.Parameters()
        params.add('EWRC', value=initial_guess[0], min=0.01, max=8.0)
        params.add('B', value=initial_guess[1], min=0., max=8.0)
        params.add('MRC', value=initial_guess[2], min=12., max=17.)
        params.add('SIGMA', value=initial_guess[3], min=0.01, max=2.0)

        mini = lmfit.Minimizer(self.log_likelihood_nataf_neg, params, fcn_args=(M, ))
        results = mini.minimize(method='nelder')

        del mini

        best_params = [results.params['EWRC'].value, results.params['B'].value, results.params['MRC'].value, results.params['SIGMA'].value]

        # Calculating the best integral
        # Preparing the scipy integrator in C
        c_vals = (ctypes.c_double * len(best_params))(*best_params)
        user_data = ctypes.cast(ctypes.pointer(c_vals), ctypes.c_void_p)
        func = LowLevelCallable(self.lib.f, user_data)

        I = integrate.quad(func, self.cmd.cm_dict['altmag'][0], self.cmd.cm_dict['altmag'][1])

        A = self.N_obs/I[0]

        N_RC = A*best_params[0]

        best_fit_params =    {'A':A, 'B':best_params[1], 'MRC':best_params[2], 'SIGMA':best_params[3], 'NRC':N_RC, 'EWRC':best_params[0]}

        # results = optimize.minimize(self.log_likelihood_nataf, (1.2,0.8,14,0.5), args=(M,), method='Nelder-Mead', bounds=((0.01, 8), (0, 4), (12, 17), (0.1, 2)))
        # ipdb.set_trace()
        return results, best_fit_params

#-----------------------------------------------------------------------
# Components to the luminosity function for the number of stars at
# magnitude K in an interval dK.
#-----------------------------------------------------------------------

def RCgaussian(M,M_RC,sigma_RC,N_RC):
    """Gaussian term for the Red Clump stars.

    Parameters
    ----------
    M : float
        Magnitude of the stars
    M_RC : float
        Mean magnitude parameter for the RC.
    sigma_RC : float
        Magnitude dispersion.
    N_RC : float
        Number of RC stars.

    Returns
    -------
    function
        Function for the gaussian term of the RC.
    """    

    return N_RC/(sqrt2pi*sigma_RC)*np.exp(-(M-M_RC)**2/(2*sigma_RC**2))

def RGBBgaussian(M,M_RC,sigma_RC,N_RC):
    """Gaussian term for the Red Giant Branch Bump stars.

    Parameters
    ----------
    M : float
        Magnitude of the stars
    M_RC : float
        Mean magnitude parameter for the RC.
    sigma_RC : float
        Magnitude dispersion.
    N_RC : float
        Number of RC stars.

    Returns
    -------
    function
        Function for the gaussian term of the RGBB.
    """    

    # Relation parameters
    N_RGBB = 0.201*N_RC
    M_RGBB = M_RC+0.737
    sigma_RGBB = sigma_RC

    return N_RGBB/(sqrt2pi*sigma_RGBB)*np.exp(-(M-M_RGBB)**2/(2*sigma_RGBB**2))

def BackgroundExp(A, B, M, M_RC):
    """The background exponential term of the luminosity function.

    Parameters
    ----------
    A : float
        Free parameter for the exponential term.
    B : float
        Free parameter for the exponential term.
    M : float
        Magnitude of the stars.
    M_RC : float
        Mean magnitude parameter for the RC.

    Returns
    -------
    function
        Function for the gaussian term of the RGBB.
    """
    return A*np.exp(B*(M-M_RC))

def BackgroundExp_NoA(B, M, M_RC):
    """The background exponential term of the luminosity function without
    the A parameter. Used in the integration technique to solve for A.

    Parameters
    ----------
    B : float
        Free parameter for the exponential term.
    M : float
        Magnitude of the stars.
    M_RC : float
        Mean magnitude parameter for the RC.

    Returns
    -------
    function
        Function for the gaussian term of the RGBB.
    """
    return np.exp(B*(M-M_RC))

def luminosityFunction(M,A,B,M_RC,sigma_RC,N_RC):
    """The luminosity function for the RC stars.

    Parameters
    ----------
    M : float
        Magnitude of the stars
    A : float
        Free parameter for the exponential term.
    B : float
        Free parameter for the exponential term.
    M_RC : float
        Mean magnitude parameter for the RC.
    sigma_RC : float
        Magnitude dispersion.
    N_RC : float
        Number of RC stars.

    Returns
    -------
    function
        Function for the luminosity function.
    """
    return RCgaussian(M,M_RC,sigma_RC,N_RC)+RGBBgaussian(M,M_RC,sigma_RC,N_RC)+BackgroundExp(A,B,M,M_RC)

def luminosityFunction_EWRC(M,EWRC,B,M_RC,sigma_RC):
    """The luminosity function for the RC stars using the equivalent width of the RC.

    Parameters
    ----------
    M : float
        Magnitude of the stars
    EWRC : float
        Equivalent width of the RC, defined as N_RC/A.
    B : float
        Free parameter for the exponential term.
    M_RC : float
        Mean magnitude parameter for the RC.
    sigma_RC : float
        Magnitude dispersion.

    Returns
    -------
    function
        Function for the luminosity function.
    """
    return RCgaussian(M,M_RC,sigma_RC,EWRC)+RGBBgaussian(M,M_RC,sigma_RC,EWRC)+BackgroundExp_NoA(B,M,M_RC)

def calcWeights(data,A,B,M_RC,sigma_RC,N_RC):
    """Calculates the weights to be used in determining the color of the RC.

    Parameters
    ----------
    A : float
        Best fit parameter for the exponential term.
    B : float
        Best fit parameter for the exponential term.
    M_RC : float
        Best fit parameter for the magnitude of the RC.
    sigma_RC : float
        Best fit parameter for the magnitude dispersion.
    N_RC : float
        Best fit parameter for the number of RC stars.

    Returns
    -------
    float
        Weight for the RC gaussian term.
    """    

    #Calculate weights using Nataf 2013 equation 6
    M = data[:,0]
    weights = RCgaussian(M,M_RC,sigma_RC,N_RC)/luminosityFunction(M,A,B,M_RC,sigma_RC,N_RC)
    return weights

def determineColor(data,weights,ColorRCinput):
    """
    Uses a brute force method to minimize the weighted variance and finds the color of the Red Clump.
    
    """
    fitMags = []
    fitHMK = []
    #Testing a min mu value of -0.5 
    MU1min = max(ColorRCinput-2.5,0)
    #MU1min = -0.5 #This is a test to see if its beuno (it's not)
    
    MU1max = ColorRCinput-.01
    MU1step = 0.01
    MU1s = np.arange(MU1min,MU1max,step=MU1step) #Makes a range of values to brute force the color over
    MU2 = ColorRCinput
    Kmags = data[:,0] #Currently not being used, pulls K-band magnitude from best fit
    HMKdiffer = data[:,1] #Pulls the best fit for (H-K) color from CMD
    VarianceMin = 1000000
    for MU1 in MU1s:
        
        #calculate differences between of H-K color and foreground disk guess and RC Color Input
        MU1diff = np.abs(HMKdiffer-MU1)
        MU2diff = np.abs(HMKdiffer-ColorRCinput)
        #find indices that are closer to colorRCinput
        #RCCinds are RCColor indices where MU2diff is less than MU1
        #actually, I'm going to carry these inds into future calcs rather than make new arrays
        #THESE ARE TRUE IF THE COLOR IS CLOSER TO THE RC
        #This returns true if a stars color is closer to the red clump than the foreground
        RCCinds = MU1diff>MU2diff
        #those colors closer to RC guess
        HMKverygood = HMKdiffer[RCCinds]
        #weights of those colors closer to RC
        #XXX still need to calc weights
        #weightsverygood = weights[RCCinds]
        #squared color differneces for RedClumb
        sum2b = HMKdiffer[RCCinds]**2
        #squared sum of color distances for those closer to the foreground disk
        sum2a = HMKdiffer[~RCCinds]**2
        #color status for fitting later, whether rejected or not
        ColorStatus = np.ones(len(HMKdiffer),dtype=bool) #list of True boolean values
        
        #loop controller
        passer = 1
        passer_counter = 0
        
        while passer==1 and passer_counter<100:
            passer = 0
            ColorRC = 0
            sumRC = 0
            combinds = np.logical_and(ColorStatus, RCCinds) #True if ColorStatus and RCCinds both true
            
            sumRC = np.sum(weights[combinds]) #Sum of weights
            sumColorRC = np.sum(weights[combinds]*HMKdiffer[combinds])#Sum of weights and color
            ColorRC = sumColorRC/sumRC #Takes the weighted average of the RC color
            
            #loop through all the color measurements
            #I think this can be done by ANDing the ColorStatus and RCCinds
            #so only include HMdiffer elements if they pass both
            V1 = np.sum(weights[combinds])
            V2 = np.sum(weights[combinds]**2)
            VARsum = np.sum(weights[combinds]*(HMKdiffer[combinds]-ColorRC)**2)
            #Now obtain the variance
            STD2 = (V1/(V1*V1 - V2))*VARsum
            #Sigma is sqrt variance
            STDcolor = np.sqrt(STD2)
            #the cutoff is 2.5 sigma
            limitcutoff = 2.5*STDcolor
            #Start sigma clipping for 2.5 sigma
            #Creates a bolean array based on whether a point is within cutoff
            #sigmainds are True for those that pass the sigma clipping (and false otherwise)
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
                
            #if passer_counter>95:
                #ipdb.set_trace()
        if passer_counter>95:
            #XXX
            #ipdb.set_trace()
            print('passer_counter limit reached')
            return [float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]
            #add logging
            
        #thus ends passer loop, back in MU1iter
        #factor of 1.096? seems familiar, can't recall why its scaling it up
        varianceB = (1.096*STDcolor)*(1.096*STDcolor);
        RealSTDcolor = 1.096*STDcolor;
        #add up the lower ninty percent of the sum2a array. If its greater than 1, calc variance 
        #otherwise just use the varaince as the calc variance
        sumsquares = np.sum(np.sort(sum2a)[:int(.9*len(sum2a))])
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
    try:
        return [FinalColorRC, RealSTDcolorOpt, STDtotal, VarianceMin, MU1opt, MU2opt]
    except:
        #XXX VarianceTotal can be nan
        #ipdb.set_trace()
        return [float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')]



if __name__ == "__main__":
    
    tstart = time.time()
    # cmd = createCMD.cmd(l=-2,b=1.5, edge_length=pram.arcmin/60.)
    # cmd = createCMD.cmd(263.585168,-29.938195, edge_length=pram.arcmin/60.)
    cmd = createCMD.cmd(265.685168,-27.238195, edge_length=pram.arcmin/60.)
    rc_fitter = RedClump(cmd, binnumber=int(.05*len(cmd.fitStarDict['altmag'])), iterations=50000)
    M_dumb, N, Nerr = rc_fitter.prepare_data_hist(rc_fitter.cmd.fitStarDict['altmag'])
    M = rc_fitter.cmd.fitStarDict['altmag']

    res, params = rc_fitter.run_minimizer(M)

    print(res)
    # ipdb.set_trace()
    # sampler, params, unc = rc_fitter.run_MCMC(M)
    tstop = time.time()
    print('Time taken: ', tstop-tstart)
    print(params)
    xval = np.linspace(12, 16, 1000)

    plt.bar(M_dumb,N, color='dimgray', width=0.7*4/len(M_dumb))
    # for theta in samples[np.random.randint(len(samples), size=100)]:
        # plt.plot(M, rc_fitter.model_MCMC_nataf(theta, M), color="r", alpha=0.1)
    best_model = rc_fitter.model_MCMC_nataf(res.params, xval)
    plt.plot(xval, best_model*4/len(M_dumb), color="k", lw=2, alpha=0.8)
    plt.xlabel('Magnitude')
    plt.ylabel(r'Number of Stars')
    plt.show()

    # ipdb.set_trace()

    

    # labels = ['EWRC','B',r'$M_{RC}$',r'$\sigma_{RC}$',r'$N_{RC}$']
    plt.clf()
    fig = plt.figure(figsize=(7, 7))
    # fig = corner.corner(samples,show_titles=True,labels=labels,plot_datapoints=True,quantiles=[0.16, 0.5, 0.84], fig=fig)
    # plt.show()
    emcee_plot = corner.corner(res.flatchain, show_titles=True, labels=res.var_names, plot_datapoints=True, quantiles=[0.16, 0.5, 0.84], fig=fig)#, truths=list(res.params.valuesdict().values()))
    plt.show()