import numpy as np
import matplotlib.pyplot as plt
import make_map as mm
import ipdb
import parameters as pram


def relUncertaintyHist(map_data,axis=10):
    ipdb.set_trace()
    header = ['RA','DEC','l','b','edgelength','N_stars','A','Aerr','B','Berr','M_RC','M_RCerr',\
        'sigma_RC','sigma_RCerr','N_RC','N_RCerr','FinalColorRC','RealSTDcolorOpt','STDtotal','VarianceMin','MU1opt','MU2opt']
    relunc = []
    for pixel in map_data:
        #relunc.append(pixel[axis+1]/pixel[axis]) #Normal uncertainty making

        #Sets any crazy high and non existing uncertainties to 1
        rel = pixel[axis+1]/pixel[axis]
        if rel > 1 or rel == 0 or np.isnan(rel) == True:
            rel = 1
        
        relunc.append(rel)

    bins = np.logspace(-4,np.log10(2),40)
    vals,ubins = np.histogram(relunc,bins=bins)
    width = 0.7*(np.log10(ubins[1])-np.log10(ubins[0]))
    xvals = (np.log10(ubins[1:])+np.log10(ubins[:-1]))/2
    ipdb.set_trace()
    print(relunc)
    print(vals,ubins)
    plt.bar(xvals,np.log10(vals),width=width)
    plt.xlabel('Fractional Uncertainty')
    plt.ylabel('Occurrences')
    plt.title('Fractional Uncertainty for '+header[axis])
    plt.show()

def generalHist(map_data,axis=14, ax=None):
    # header = ['RA','DEC','l','b','edgelength','N_stars','A','Aerr','B','Berr','M_RC','M_RCerr',\
        # 'sigma_RC','sigma_RCerr','N_RC','N_RCerr','FinalColorRC','RealSTDcolorOpt','STDtotal','VarianceMin','MU1opt','MU2opt']
    header = ['RA','DEC','l','b','A','B','M_RC',\
        'sigma_RC','N_RC','FinalColorRC','RealSTDcolorOpt','STDtotal','VarianceMin','MU1opt','MU2opt']

    values = [] #list to isolate specific axis
    for pixel in map_data:
        values.append(pixel[axis])
    # ipdb.set_trace()
    # bins = np.linspace(np.nanmin(values),np.nanmax(values),31) 
    bins = np.linspace(0,1,31) 
    vals,ubins = np.histogram(values,bins=bins)
    width = 0.7*(ubins[1]-ubins[0])
    xvals = (ubins[1:]+ubins[:-1])/2
    

    ax.bar(xvals,vals,width=width)
    ax.vlines(np.nanmedian(values),0,np.nanmax(vals),color='r',linestyle='--',label='Median='+str(round(np.nanmedian(values),3)))
    ax.vlines(np.quantile(values,0.16),0,np.nanmax(vals),color='g',linestyle='--',label='16th Percentile='+str(round(np.quantile(values,0.16),3)))
    ax.vlines(np.quantile(values,0.84),0,np.nanmax(vals),color='g',linestyle='--',label='84th Percentile='+str(round(np.quantile(values,0.84),3)))
    ax.legend()
    ax.set_xlabel(header[axis])
    ax.set_ylabel('Occurrences')
    ax.set_title('Histogram of total '+header[axis])
    # plt.savefig('figs/Total NRC Histogram.png')
    # plt.show()

if __name__ == "__main__":
    # data = mm.read_map('maps/mcmc_map_prop.map')
    fig, axs = plt.subplots(1,1)
    # axs = axs.flatten()
    # data = mm.read_map('maps/map_'+pram.phot+'_'+str(pram.year)+'_'+str(pram.arcmin)+'.map')
    data = mm.read_map('maps/mcmc_map_1.5.map')
    # for i, ax in enumerate(axs):

    generalHist(data, axis=7, ax=axs)
    plt.show()
