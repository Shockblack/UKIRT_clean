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
from scipy.special import erf
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
plt.style.use('az-paper-twocol')

from multiprocessing import Pool

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

# B = 0.43
sqrt2pi = np.sqrt(2*np.pi)

class RedClump():
    def __init__(self, cmd, nwalkers=50, iterations=1000, burnin=100):
        self.cmd = cmd
        # self.binnumber = int(.05*len(cmd.fitStarDict['altmag']))
        self.binnumber = 40
        self.nwalkers = nwalkers
        self.iterations = iterations
        self.burnin = burnin

        self.N_obs = len(cmd.fitStarDict['altmag'])

        #-------------DEPRECATED----------------#
        # Importing the compiled C function to integrate
        # self.lib = ctypes.CDLL(os.path.abspath('integratelib.so'))
        # self.lib.f.restype = ctypes.c_double
        # self.lib.f.argtypes = (ctypes.c_int, ctypes.POINTER(ctypes.c_double), ctypes.c_void_p)
        #---------------------------------------#


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
        EWRC, M_RC, sigma_RC, B = theta#['EWRC'].value, theta['B'].value, theta['MRC'].value, theta['SIGMA'].value
        
        A = self.A

        N_RC = EWRC*A

        # Move log space to linear

        return luminosityFunction(M,A,B,M_RC,sigma_RC,N_RC)
    
    def model_lmfit(self, theta, M, method=None):
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
        try:
            EWRC, M_RC, sigma_RC, B = theta['EWRC'].value, theta['MRC'].value, theta['SIGMA'].value, theta['B'].value
        except:
            pass
        try:
            EWRC, M_RC, sigma_RC , B= theta['EWRC'], theta['MRC'], theta['SIGMA'], theta['B']
        except:
            pass
        try:
            EWRC, M_RC, sigma_RC, B = theta[0], theta[1], theta[2], theta[3]
        except:
            pass
        # Integrate
        I = self.integrator(EWRC, B, M_RC, sigma_RC)
        A = self.N_obs/I

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
        EWRC, M_RC, sigma_RC, B = theta
        if (0. < EWRC < 10. and 12. < M_RC < 18. and 0. < sigma_RC < 1. and 0. < B < 1.):
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

        EWRC, M_RC, sigma_RC, B = theta

        self.I = self.integrator(EWRC, B, M_RC, sigma_RC)
        self.A = self.N_obs/self.I
        model = self.model_MCMC_nataf(theta, M)

        return np.sum(np.log(model))-self.N_obs

    def log_likelihood_nataf_neg(self, theta, M):
        # import ipdb; ipdb.set_trace()
        
        EWRC, M_RC, sigma_RC, B = theta['EWRC'].value, theta['MRC'].value, theta['SIGMA'].value, theta['B'].value

        # ipdb.set_trace()
        # c_vals = (ctypes.c_double * len(py_vals))(*py_vals)
        # user_data = ctypes.cast(ctypes.pointer(c_vals), ctypes.c_void_p)
        # func = LowLevelCallable(self.lib.f, user_data)

        # I = integrate.quad(func, self.cmd.cm_dict['altmag'][0], self.cmd.cm_dict['altmag'][1])
        # if I[0]==0:
        #     ipdb.set_trace()
        # self.I = I[0]
        self.I = self.integrator(EWRC, B, M_RC, sigma_RC)
        self.A = self.N_obs/self.I
        model = self.model_lmfit(theta, M)

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
        # bins = np.linspace(self.cmd.cm_dict['altmag'][0], self.cmd.cm_dict['altmag'][1],self.binnumber)
        # Calculate number of bins from data
        num_bins = int(np.sqrt(len(M)))+1
        bins = np.linspace(min(M), max(M), num_bins)

        Nstars, hist_bins = np.histogram(M,bins)#'fd')

        mags = (hist_bins[1:]+hist_bins[:-1])/2
        
        return mags, Nstars

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
            initial_guess = np.array([2.,14,0.5,0.43])
        # if initial_guess[3] > 0.49:
        #     initial_guess[3] = 0.45
        
        self.initial_guess = initial_guess
        # p0 = [np.array(initial_guess) + 0.05*np.array(initial_guess)*np.random.randn(len(initial_guess)) for i in range(self.nwalkers)]
        p0 = [np.array(initial_guess) + 10**(-5)*np.array(initial_guess)*np.random.randn(len(initial_guess)) for i in range(self.nwalkers)]

        # with Pool(4) as pool:

        sampler = emcee.EnsembleSampler(self.nwalkers, len(initial_guess), self.log_probability, args=(M, ))#, pool=pool)

        # Run the burn-in
        p0, _, _ = sampler.run_mcmc(p0, self.burnin, progress=False)
        sampler.reset()

        autocorr = np.empty(self.iterations)

        # This will be useful to testing convergence
        old_tau = np.inf
        i = 0
        # Now we'll sample for up to max_n steps
        for sample in sampler.sample(p0, iterations=self.iterations, progress=True):
            # Only check convergence every 200 steps
            if sampler.iteration % 200:
                continue
            
            # Compute the autocorrelation time so far
            # Using tol=0 means that we'll always get an estimate even
            # if it isn't trustworthy
            tau = sampler.get_autocorr_time(tol=0)
            
            
            autocorr[i] = np.mean(tau)
            i += 1
            
            # Check convergence
            converged = np.all(tau * 100 < sampler.iteration)
            converged &= np.all(np.abs(old_tau - tau) / tau < 0.01)
            if converged:
                break
            old_tau = tau
            best_ind = np.argmax(sampler.flatlnprobability)
            samples = sampler.flatchain
            best_params = samples[best_ind]
            
        # Run the production
        # pos, prob, state =  sampler.run_mcmc(p0, self.iterations, progress=True)

        # Get the best fit parameters
        best_ind = np.argmax(sampler.flatlnprobability)
        samples = sampler.flatchain
        best_params = samples[best_ind]
        B = best_params[3]
        # Calculate the statistical uncertainties
        perc = np.percentile(samples, [16, 50, 84], axis=0)
        unc = np.diff(perc, axis=0).T
        unc = np.sqrt(unc[:,0]**2 + unc[:,1]**2)

        # Calculating the best integral
        # Preparing the scipy integrator in C
        # c_vals = (ctypes.c_double * len(best_params))(*best_params)
        # user_data = ctypes.cast(ctypes.pointer(c_vals), ctypes.c_void_p)
        # func = LowLevelCallable(self.lib.f, user_data)

        # I = integrate.quad(func, self.cmd.cm_dict['altmag'][0], self.cmd.cm_dict['altmag'][1])
        self.I = self.integrator(best_params[0], B, best_params[1], best_params[2])
        A = self.N_obs/self.I#I[0]
        self.A = A
        self.A_mc = A
        # Propogating the error
        # A_err = A*np.sqrt( (np.sqrt(self.N_obs)/self.N_obs)**2 + (I[1]/I[0])**2 )
        N_RC = A*best_params[0]
        # N_RC_err = N_RC*np.sqrt((A_err/A)**2 + (unc[0]/best_params[0])**2)
        A_err = 0.
        N_RC_err = 0.

        params =    {'A':A, 'A_err':A_err, 'B':B, 'MRC':best_params[1], 'MRC_err':unc[1], \
                    'SIGMA':best_params[2], 'SIGMA_err':unc[2], 'NRC':N_RC, 'NRC_err':N_RC_err, 'EWRC':best_params[0], 'EWRC_err':unc[0]}

        self.mcmc_params = {'EWRC':best_params[0], 'B':B, 'MRC':best_params[1], 'SIGMA':best_params[2]}
        self.samples = sampler.flatchain
        return sampler, params, unc


    def integrator(self, EWRC, B, M_RC, sigma_RC):
        """
        Integrates the luminosity function to find the normalization constant A
        by solving the definite integral and plugging in the bounds of the fit data.
        """
        min_K = self.cmd.cm_dict['altmag'][0]
        max_K = self.cmd.cm_dict['altmag'][1]

        term1 = BackgroundExp_NoA(B, max_K, M_RC)/B - BackgroundExp_NoA(B, min_K, M_RC)/B
        term2 = EWRC/2 * (erf((max_K-M_RC)/(np.sqrt(2)*sigma_RC)) - erf((min_K-M_RC)/(np.sqrt(2)*sigma_RC)))
        EWRGBB = 0.201*EWRC
        M_RGBB = M_RC+0.737
        term3 = EWRGBB/2 * (erf((max_K-M_RGBB)/(np.sqrt(2)*sigma_RC)) - erf((min_K-M_RGBB)/(np.sqrt(2)*sigma_RC)))

        self.integral = term1 + term2 + term3
        return self.integral


    def run_minimizer(self, M, initial_guess=None):

        self.value_tracker = []

        if type(initial_guess) == type(None):
            initial_guess = np.array([2.,14,0.5])


        params = lmfit.Parameters()
        params.add('EWRC', value=initial_guess[0], min=0.01, max=20.0)
        params.add('MRC', value=initial_guess[1], min=12., max=18.)
        params.add('SIGMA', value=initial_guess[2], min=0.01, max=1.)
        params.add('B', value=0.43, min=0.01, max=1.0)
        
        mini = lmfit.Minimizer(self.log_likelihood_nataf_neg, params, fcn_args=(M, ))
        results = mini.minimize(method='nelder')
        self.results = results
        del mini

        best_params = [results.params['EWRC'].value, results.params['MRC'].value, results.params['SIGMA'].value]
        B = results.params['B'].value
        # Calculating the best integral
        # Preparing the scipy integrator in C
        # c_vals = (ctypes.c_double * len(best_params))(*best_params)
        # user_data = ctypes.cast(ctypes.pointer(c_vals), ctypes.c_void_p)
        # func = LowLevelCallable(self.lib.f, user_data)

        # I = integrate.quad(func, self.cmd.cm_dict['altmag'][0], self.cmd.cm_dict['altmag'][1])

        # A = self.N_obs/I[0]
        I = self.integrator(best_params[0], B, best_params[1], best_params[2])
        self.A = self.N_obs/I
        N_RC = self.A*best_params[0]
        
        self.I_nm = I
        self.A_nm = self.N_obs/I
        self.N_RC_nm = N_RC

        best_fit_params =    {'A':self.A, 'B':B, 'MRC':best_params[1], 'SIGMA':best_params[2], 'NRC':N_RC, 'EWRC':best_params[0]}

        self.best_fit_params = best_fit_params

        return results, best_fit_params

    def plot(self, ax=None, fig=None, show=False):
        # Make the figure
        if type(ax) == type(None):
            fig, ax = plt.subplots(1,1)
        M_dumb, N = self.prepare_data_hist(self.cmd.fitStarDict['altmag'])
        
        star_max, star_min = max(self.cmd.fitStarDict['altmag']), min(self.cmd.fitStarDict['altmag'])

        xval = np.linspace(star_min, star_max, 1000)

        x_range = abs(star_max - star_min)

        ax.bar(M_dumb, N, color='dimgray', width=0.7*x_range/len(M_dumb))

        initial_model = self.model_lmfit(self.results.params, xval, method='nm')
        best_model = self.model_lmfit(self.mcmc_params, xval, method='mc')

        # Add text in the top right corner indicating the location of the field
        ax.text(0.95, 0.95, f'(l,b)={self.cmd.l, self.cmd.b}', fontsize=14, va='top', ha='right', color='k', weight='bold', transform=ax.transAxes)


        model_list = []

        # Loop over samples and calculate models
        for theta in self.samples:
            model_list.append(self.model_lmfit(theta, xval))

        # Get the quantiles for the model
        quantiles = np.percentile(model_list, [2.5, 16, 50, 84, 97.5], axis=0)
        
        # Plot the quantiles using fill_between
        ax.fill_between(xval, quantiles[0]*x_range/len(M_dumb), quantiles[4]*x_range/len(M_dumb), color=pram.red, alpha=0.2)
        ax.fill_between(xval, quantiles[1]*x_range/len(M_dumb), quantiles[3]*x_range/len(M_dumb), color=pram.red, alpha=0.5)

        # ax.plot(xval, quantiles[2]*x_range/len(M_dumb), color="r", lw=2, alpha=0.8, label='MCMC')
        ax.plot(xval, best_model*x_range/len(M_dumb), color='r', lw=4, alpha=1.0, label='MCMC')
        ax.plot(xval, initial_model*x_range/len(M_dumb), color="k", lw=2, alpha=0.8, label='Nelder-Mead', ls='--')
        ax.set_xlim(star_min, star_max)
        
        ax.set_xlabel(r'$K$-Band Magnitude', fontsize=16)
        ax.set_ylabel(r'Number of Stars', fontsize=16)
        ax.legend(fontsize=14, loc='upper left', handlelength=1.2)

        ax.tick_params(axis='both',which='both',direction='in',top=True,right=True, labelsize=16)

        if show:
            plt.show()

    def plotall(self, ax=None, fig=None, show=False):

        # Make the figure
        if type(ax) == type(None):
            fig, ax = plt.subplots(1,1)
        M_dumb, N = self.prepare_data_hist(self.cmd.filterStarDict['altmag'])
        star_max, star_min = max(self.cmd.filterStarDict['altmag']), min(self.cmd.filterStarDict['altmag'])

        xval = np.linspace(star_min, star_max, 1000)

        x_range = abs(star_max - star_min)

        ax.bar(M_dumb, N, color='dimgray', width=0.7*x_range/len(M_dumb))

        initial_model = self.model_lmfit(self.results.params, xval, method='nm')
        best_model = self.model_lmfit(self.mcmc_params, xval, method='mc')

        # Add text in the top right corner indicating the location of the field
        ax.text(0.99*star_max, max(N), f'(l,b)={round(self.cmd.l, 3), round(self.cmd.b, 3)}', fontsize=12, va='top', ha='right', color='k', weight='bold')

        for theta in self.samples[np.random.randint(len(self.samples), size=100)]:
            ax.plot(xval, self.model_lmfit(theta, xval)*x_range/len(M_dumb), color="r", lw=1, alpha=0.1)

        ax.plot(xval, best_model*x_range/len(M_dumb), color="r", lw=2, alpha=0.8, label='MCMC')
        ax.plot(xval, initial_model*x_range/len(M_dumb), color="k", lw=2, alpha=0.8, label='Nelder-Mead', ls='--')
        ax.set_xlim(star_min, star_max)
        
        ax.set_xlabel(r'$K$-Band Magnitude', fontsize=16)
        ax.set_ylabel(r'Number of Stars', fontsize=16)
        ax.legend(fontsize=14, loc='upper left', handlelength=1.5)

        ax.tick_params(axis='both',which='both',direction='in',top=True,right=True, labelsize=16)

        if show:
            plt.show()

    def plot_full(self, show=False):
        # Need to use new matplotlib feature, subfigures
        fig = plt.figure(figsize=(12,4))
        subfigs = fig.subfigures(1, 2, wspace=0, width_ratios=[2, 1])

        # Add the model fit
        ax1 = subfigs[0].subplots(1,1)

        subfigs[0].subplots_adjust(wspace=0.)
        subfigs[1].subplots_adjust(wspace=0.)

        self.plot(ax=ax1, fig=subfigs[0], show=False)

        ax1.set_xlim(12, 16.5)

        # Add the corner plot
        labels = [r'$EW_{RC}$', r'$K_{RC}$', r'$\sigma_{RC}$', r'$B$']
        emcee_plot = corner.corner(rc_fitter.samples, show_titles=True, labels=labels, plot_datapoints=True, quantiles=[0.16, 0.5, 0.84], truths=list(res.params.valuesdict().values()), fig=subfigs[1])


        from matplotlib.ticker import FormatStrFormatter
        for i, ax in enumerate(emcee_plot.get_axes()):
            ax.title.set_fontsize(10)
            ax.set_xlabel(ax.get_xlabel(), fontsize=10)
            ax.set_ylabel(ax.get_ylabel(), fontsize=10)
            for label in ax.get_xticklabels() + ax.get_yticklabels():
                label.set_fontsize(10)
            ax.tick_params(pad=1.5)
            ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
            ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))

            if i < 6:
                ax.xaxis.set_ticklabels([])
            if i % 3 != 0:
                ax.yaxis.set_ticklabels([])
                
        if show:
            plt.show()

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

def AGBBgaussian(M,M_RC,sigma_RC,N_RC):

    N_AGBB = 0.028*N_RC
    M_AGBB = M_RC-1.07
    sigma_AGBB = sigma_RC

    return N_AGBB/(sqrt2pi*sigma_AGBB)*np.exp(-(M-M_AGBB)**2/(2*sigma_AGBB**2))

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
    return RCgaussian(M,M_RC,sigma_RC,EWRC)+RGBBgaussian(M,M_RC,sigma_RC,EWRC)+BackgroundExp_NoA(B,M,M_RC)+AGBBgaussian(M,M_RC,sigma_RC,EWRC)

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
    l = 0.5
    b = -0.25
    # b=-2.0
    cmd = createCMD.cmd(l=l,b=b, edge_length=pram.arcmin/60.)
    # cmd = createCMD.cmd(263.585168,-29.938195, edge_length=pram.arcmin/60.)
    # cmd = createCMD.cmd(265.685168,-27.238195, edge_length=pram.arcmin/60.)
    rc_fitter = RedClump(cmd, iterations=100000)
    M_dumb, N = rc_fitter.prepare_data_hist(rc_fitter.cmd.fitStarDict['altmag'])
    M = rc_fitter.cmd.fitStarDict['altmag']

    res, params = rc_fitter.run_minimizer(M)

    init_params = np.array([params['EWRC'], params['MRC'], params['SIGMA'], params['B']])

    # ipdb.set_trace()
    sampler, params, unc = rc_fitter.run_MCMC(M, initial_guess=init_params)
    tstop = time.time()
    print('Time taken: ', tstop-tstart)


    rc_fitter.plot_full(show=True)

    # plt.savefig('paperfigs/fit_with_corner_l'+str(l)+'_b'+str(b)+'.pdf', bbox_inches='tight')