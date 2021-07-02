def setOptions():

  options={
    # path to the lightcurve data
    #'dataDir':'/data/microlensing',  
    #'dataDir':'//Users/johnson.7080/work/jplsip/extinctionMaps',  #SAMSON
    #Aiden's directory below
    'dataDir':'/mnt/a/documents/files/surp/data', 
    # PSF data is too big to stote locally; define an alternate data directory
    'dataDirPSF':'/Volumes/Data/',  
    
    #'photomType':'CASU',     # choose CASU or PSF data reduction
    'photomType':'PSF',      # choose CASU or PSF data reduction
    'year':'2016',           # select which year of observation
    #'year':'2015',           # select which year of observation
    'year':'2017',           # select which year of observation
    #'year':'2018',           # select which year of observation
    'field':'s2_1',          # sometimes you can select an individual field here
    #'field':'11',          # sometimes you can select an individual field here
    'ccd':'3',               # sometimes you can select an individual CCD here

    'overlapData':0,         # use the full dataset for each lightcurve?
    'overlapData':1,         # use the full dataset for each lightcurve?
    'selectionData':0,       # use the original dChi2 for selection or overlap?
    #'selectionData':1,       # use the original dChi2 for selection or overlap?

    #'deltaChiThreshold':1000,  # used as threshold in mcmcFit,slopeFit,sinFit
    'deltaChiThreshold':500,  # used as threshold in mcmcFit,slopeFit,sinFit
    #'deltaChiThreshold':300,  # used as threshold in mcmcFit,slopeFit,sinFit
    #'deltaChiThreshold':100,  # used as threshold in mcmcFit,slopeFit,sinFit
    #'deltaChiThreshold':101,  # used as threshold in mcmcFit,slopeFit,sinFit
    #'deltaChiThreshold':50,  # used as threshold in mcmcFit,slopeFit,sinFit
    #'deltaChiThreshold':313,  # used as threshold in mcmcFit,slopeFit,sinFit
    #'deltaChiThreshold':666,  # used as threshold in mcmcFit,slopeFit,sinFit
    
    'gouldStyle':1,           # use the faster Gould/Kim 2017 method for the first fits

    'startFromPreviousFit':0, # use previous fit as starting parameters for MCMC

    'JDadjust':2450000.,      # subtract this from HJD times

    'justUKIRT':1,            # just fit the UKIRT data
    'includeSpitzer':0,       # also consider Spitzer data
    'includeMOA':0,           # also consider MOA data

    'includeBlending':1,      # always do this (otherwise f_b=1)
    'properBlending':1,       # analytic solution for F_s and F_b

# remove 5-sigma outliers if more than 3 t_E away from t_0
    'baselineCut':3.,
    'sigmaClip':5.,

#    'nTimeCuts':31,
#    'nTimeCuts':5,
    'nTimeCuts':1,
    }

# parameters for setting a grid-based fit (Gould-style)
  gridOptions={
# range and spacing for the t_eff grid:
    'teffmin':1, 'teffmax':45, 'teffFac':4./3.,
# t_0 spacing: in units of t_eff
    't0stepFac':1./3.,
# t_0 range: number of t_eff's outside data range
    't0rangeFac':1./2.
    }

#____________________MCMC Options____________________
  hires=1
  hires=0
  if hires:
    nsteps0=1000
    nsteps=500
  else:
    nsteps0=100
    nsteps=50

# this is better for fitting fake data,
#  but really must match exactly the real data analysis
#  (hold on - they are different, as far as blending fitting)
#    nsteps0=300
#    nsteps=100

# 100/50 is not fully converged on a small-u0 fakedata run
# but maybe this is o.k.; still has a very good delta-chi2

# for testing:
#    nsteps0=1
#    nsteps=500

  mcmcOptions={
    'nwalkers':20,          # number of MCMC walkers
    'nsteps0':nsteps0,      # number of MCMC steps for initial burn-in
    'nsteps':nsteps         # number of MCMC steps (after burn-in)
    }

  plotOptions={
    'lightcurvePlot':1,       # plot lightcurves and model fits?
    'trianglePlot':0,         # plot the correlations between variables?
    'markertrendPlot':0       # plot the chain evolution?
    }

  return options,gridOptions,mcmcOptions,plotOptions
#________________________________________________________________________________________


