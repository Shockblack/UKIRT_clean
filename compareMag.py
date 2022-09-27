#------------------------------------------------------------------------------
# Filename: compareMap.py
# 
# Compare the our maps and that of Gonzalez et al. 2012 and 2018.
#
# Programmer: Aiden Zelakiewicz
#
# Revision History:
#   02-Sep-2022 : File Created
#------------------------------------------------------------------------------

# Importing necessary packages and files
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as pat
from matplotlib import cm

import seaborn as sns

import plotly.express as px

from astropy import coordinates as coord
from astropy import units as u 

#------------------------------------------------------------------------------

def read_map(filename):
    """Importing the map data given a csv filename

    Parameters
    ----------
    filename : str
        The name of the csv file with map data.

    Returns
    -------
    map_data : numpy.ndarray
        The map data in the form of a numpy array.
    """
    map_data = np.loadtxt(filename,delimiter=',')
    return map_data

#------------------------------------------------------------------------------

def plot_point_map(map_data,axis=6):
    fig, ax = plt.subplots()
    
    # Building the color bar
    good_data = map_data[:,axis][~np.isnan(map_data[:,axis])]
    cmap_max = np.max(good_data)
    cmap_min = np.min(good_data)
    
    # Adding cmap
    norm = matplotlib.colors.Normalize(vmin=cmap_min,vmax=cmap_max)
    cmap = matplotlib.cm.get_cmap('cividis')

    sc = ax.scatter(map_data[:,0],map_data[:,1],marker='.',c=map_data[:,axis],cmap=cmap,vmin=cmap_min,vmax=cmap_max)
    fig.colorbar(sc,ax=ax)

    ax.set_aspect('equal')

#------------------------------------------------------------------------------

def calc_angle(ra1,dec1,ra2,dec2):
    c1 = coord.SkyCoord(ra=ra1*u.degree,dec=dec1*u.degree,frame='icrs')
    c2 = coord.SkyCoord(ra=ra2*u.degree,dec=dec2*u.degree,frame='icrs')
    angle = np.degrees(np.arctan2(c2.galactic.b.degree-c1.galactic.b.degree,c2.galactic.l.degree-c1.galactic.l.degree))
    return angle

#------------------------------------------------------------------------------

def do_nothing(array):
    return array 

#------------------------------------------------------------------------------

def rel_diff(array, cap=True, rel=False, adjust_filter=True, ccm=False):
    gonzData = np.loadtxt('../data/gonz/gonzData_2018.dat',delimiter=',')
    M_RC = 13.1

    if ccm:
        ak = array*1.5 # Reddening law from Cardelli
    else:
        ak = array-M_RC

    diff = []
    
    if adjust_filter:
        gonzData[:,4] = (gonzData[:,4] + 0.689*0.012)/1.069

    gonz_nonull=gonzData[:,4][~np.isnan(array)]
    ak=ak[~np.isnan(array)]

    for i in range(len(gonz_nonull)):
        # Finding the difference between the AK values
        if cap:
            ukirt_ak = np.minimum(ak[i], 1.5)
            gonz_ak = np.minimum(gonz_nonull[i], 1.5)
        else:
            ukirt_ak = ak[i]
            gonz_ak = gonz_nonull[i]
        
        # Calculates relative difference if rel=True, else calculates absolute difference
        if cap and ukirt_ak==1.5 and gonz_ak==1.5:
            # Passes if the AK values are both capped at 1.5...
            # probably not the best way to do this, but it works.
            pass
        else:
            if rel:
                diff.append(abs(ukirt_ak - gonz_ak)/ukirt_ak)
            else:
                diff.append(ukirt_ak - gonz_ak)

    diff = np.array(diff)

    return diff

#------------------------------------------------------------------------------

cmap_label = r'$A(K)_{UKIRT} - A(K_s)_{VVV}$'

def plot_grid_map(map_data,lb=True,func=do_nothing,axis=7,cb_label=cmap_label,path='figs/',figname='map.pdf',useangle=True):
    fig, ax = plt.subplots()
    
    rel = False

    #build color bar
    good_data = func(map_data[:,axis][~np.isnan(map_data[:,axis])], rel=rel)
    cmap_max = np.max(good_data)
    cmap_min = np.min(good_data)
    
    cmap_max = np.percentile(good_data,95) 
    cmap_min = np.percentile(good_data,5) 


    norm = matplotlib.colors.Normalize(vmin=cmap_min,vmax=cmap_max)
    cmap = cm.get_cmap('cividis_r')
    #sc = ax.scatter(map_data[:,0],map_data[:,1],facecolor=cmap(norm(map_data[:,axis])))
    sc = ax.scatter(map_data[:,2]+1e3,map_data[:,3],marker='.',c=func(map_data[:,axis],rel=rel),cmap=cmap,vmin=cmap_min,vmax=cmap_max)
    cbar = fig.colorbar(sc,ax=ax)
    cbar.set_label(cb_label)
    
    #Check location of grids where line discontinuity occurs
    #We first find location then try to find avg mag difference
    

        
    if lb:
        #prep for l,b
        pixel = map_data[0]
        angle = calc_angle(pixel[0]-pixel[4]/2.,pixel[1]-pixel[4]/2.,pixel[0]-pixel[4]/2.,pixel[1]+pixel[4]/2.)
        

        if pixel[2] < 358 and pixel[2]>160:
            print(pixel[0],pixel[1],pixel[2],pixel[3])

        for i in range(len(map_data)):
            pixel = map_data[i]

            if np.abs(pixel[3]) <1 and pixel[axis]<13.5:
                print(pixel[0],',',pixel[1],pixel[2],pixel[3])

            if useangle == True:
                if pixel[2]>160:
                    rect = pat.Rectangle((pixel[2]-360,pixel[3]),pixel[4],pixel[4],angle=angle,facecolor=cmap(norm(good_data[i])))
                else:
                    rect = pat.Rectangle((pixel[2],pixel[3]),pixel[4],pixel[4],angle=angle,facecolor=cmap(norm(good_data[i])))
            else:
                if pixel[2]>160:
                    rect = pat.Rectangle((pixel[2]-360,pixel[3]),pixel[4],pixel[4],facecolor=cmap(norm(good_data[i])))
                else:
                    rect = pat.Rectangle((pixel[2],pixel[3]),pixel[4],pixel[4],facecolor=cmap(norm(good_data[i])))
                
            ax.add_patch(rect)


        ax.set_xlim(2.8,-2.2)
        ax.set_ylim(-2.8,2.1)
        ax.set_xlabel('$l$')
        ax.set_ylabel('$b$')
    
    ax.set_aspect('equal')
    # plt.savefig(path+figname)
    plt.show()

#------------------------------------------------------------------------------

def plot_differences_hist(map_data,func=do_nothing,axis=7,path='figs/',figname='mag_differences.pdf',ccm=False):

    if ccm:
        diff_data = func(map_data[:,axis], rel=False, cap=False, adjust_filter=True, ccm=True)
        # Rejecting outliers outside of 2.5 sigma
        # diff_data = reject_outliers(diff_data, m=2.5)
        bins = np.linspace(-0.475,0.575, 20)
        hist, bin_edge = np.histogram(diff_data, bins=bins)

    else:
        diff_data = func(map_data[:,axis], rel=False, cap=True, adjust_filter=True)
        # Rejecting outliers outside of 2.5 sigma
        # diff_data = reject_outliers(diff_data, m=2.5)
        bins = np.linspace(-0.725,0.225, 20)
        hist, bin_edge = np.histogram(diff_data, bins=bins)

    # Doing some statistics :D
    std = np.std(diff_data)
    mean = np.mean(diff_data)
    median = np.median(diff_data)
    median_std = median/np.sqrt(len(diff_data))

    z1_left = np.percentile(diff_data, 34)
    z1_right = np.percentile(diff_data, 68)
    z2_left = np.percentile(diff_data, 2.5)
    z2_right = np.percentile(diff_data, 97.5)
    
    bin_centers = 0.5 * (bin_edge[:-1] + bin_edge[1:])
    width = 0.9*(bin_edge[1]-bin_edge[0])

    plt.style.use('bmh')
    fig, ax = plt.subplots(figsize=(8,5))
    ax.set_axisbelow(True)
    
    ax.bar(bin_centers, hist, align='center', width=width, color='cornflowerblue', label='_nolegend_')
    # fig = px.bar(x=bin_edges, y=hist)

    # Adding statistic intervals
    plt.axvline(x=median, label='Median = {}±{}'.format(round(median, 5), round(median_std, 5)), ls='-.',color='black')
    plt.axvline(x=z1_left, label='34th Percentile = {}'.format(round(z1_left, 3)), ls='--', color='black')
    plt.axvline(x=z1_right, label='68th Percentile = {}'.format(round(z1_right, 3)), ls='--', color='red')
    plt.axvline(x=z2_left, label='2.5th Percentile = {}'.format(round(z2_left, 3)), ls=':', color='black')
    plt.axvline(x=z2_right, label='97.5th Percentile = {}'.format(round(z2_right, 3)), ls=':', color='red')
    
    plt.legend(loc='upper right', prop={'size': 8})

    plt.xlabel(r'$A(K)_{UKIRT} - A(K)_{VVV}$')
    plt.ylabel('Number of Pixels')

    # fig.update_layout(
        # bargap=0.05, xaxis_title=r'$\LARGE{A(K)_{UKIRT} - A(K)_{VVV}}$', 
        # yaxis_title='Number of Pixels',
        # font=dict(
            # size=24
            # )
        # )
    
    # fig.write_image(path+figname)
    plt.show()

#------------------------------------------------------------------------------

def plot_mags(map_data, func=do_nothing, axis=10, path='figs/',figname='mag_comparison_fit.pdf'):
    gonzData = np.loadtxt('../data/gonz/gonzDataScaled.dat',delimiter=',')
    
    gonzData[:,4] = (gonzData[:,4] + 0.689*0.012)/1.069
    M_RC = 13.1
    ak = map_data[:,10]-M_RC

    # ukirt_ak = np.minimum(ak, 1.5)
    # gonz_ak = np.minimum(gonzData[:,4], 1.5)

    gonz_cap = gonz_ak!=1.5
    ukirt_cap = ukirt_ak!=1.5

    non_capped = np.logical_or(gonz_cap, ukirt_cap)

    ukirt_ak = ukirt_ak[non_capped]
    gonz_ak = gonz_ak[non_capped]

    plt.style.use('ggplot')
    # ipdb.set_trace()
    plt.scatter(ukirt_ak,gonz_ak,s=0.1,color='black')
    plt.plot([-0.3,1.5],[-0.3,1.5],color='red', linestyle='--')
    plt.xlabel(r'$A(K)_{UKIRT}$')
    plt.ylabel(r'$A(K)_{VVV}$')
    plt.xlim(-0.1,1.6)
    plt.ylim(-0.1,1.6)
    plt.show()

#------------------------------------------------------------------------------

def plot_mag_card(map_data, path='figs/',figname='mag_comparison_ccm.pdf', error_sample=False, cap=False, cap_value=1.5):

    gonzData = np.loadtxt('../data/gonz/gonzDataScaled.dat',delimiter=',')
    
    # Adjusting the 2MASS Ks magnitude to match the UKIRT K
    # extinction and error values
    gonzData[:,4] = (gonzData[:,4] + 0.689*0.012)/1.069
    gonzData[:,5] = gonzData[:,5]/1.069

    # Adjusting UKIRT H-K color to match 2MASS H-K_s
    # map_data[:,16] = map_data[:,16]*1.062 + 0.017
    # map_data[:,17] = map_data[:,17]*1.062

    # Calculating the extinction using Cardelli 1989 law
    ak = map_data[:,16] * 1.5
    ak_error = map_data[:,17] * 1.5

    plt.style.use('ggplot')

    if cap:
        ukirt_ak = np.minimum(ak, cap_value)
        gonz_ak = np.minimum(gonzData[:,4], cap_value)

        gonz_cap = gonz_ak < (cap_value - np.finfo(np.float32).eps)
        ukirt_cap = ukirt_ak < (cap_value- np.finfo(np.float32).eps)

        non_capped = np.logical_or(gonz_cap, ukirt_cap)

        ukirt_ak = ukirt_ak[non_capped]
        gonz_ak = gonz_ak[non_capped]
        ak_error = ak_error[non_capped]
        gonz_error = gonzData[:,5][non_capped]
        plt.xlim(0.2,cap_value+0.1)
        plt.ylim(0.2,cap_value+0.1)

        plt.plot([-0.1,cap_value],[-0.1,cap_value],color='red', linestyle='--')
    else:
        ukirt_ak = ak
        gonz_ak = gonzData[:,4]
        gonz_error = gonzData[:,5]

        plt.plot([-0.1,4],[-0.1,4],color='red', linestyle='--')

    ak_df = pd.DataFrame({'UKIRT':ukirt_ak, 'VVV':gonz_ak, 'UKIRT_error':ak_error, 'VVV_error':gonz_error})
    ak_df_sample = ak_df.sample(n=100)

    
    # ipdb.set_trace()
    if error_sample:
        plt.errorbar(ak_df_sample['UKIRT'],ak_df_sample['VVV'],yerr=ak_df_sample['VVV_error'],xerr=ak_df_sample['UKIRT_error'],fmt='.',color='black',markersize=1)
    else:
        plt.scatter(ukirt_ak,gonz_ak,s=0.1,color='black')
    
    plt.xlabel(r'$A(K)_{UKIRT}$')
    plt.ylabel(r'$A(K)_{VVV}$')

    plt.show()

#------------------------------------------------------------------------------

def plot_contour(map_data, path='figs/',figname='mag_contour_ccm.pdf', datafile='../data/gonz/gonzDataScaled.dat', ccm=True):
    """Plots the contour map of UKIRT vs VVV extinction values using the Cardelli 1989 law

    Parameters
    ----------
    map_data : np.array
        The map data of the UKIRT data.
    path : str, optional
        The location to save the figure, by default 'figs/'
    figname : str, optional
        The name of the figure including type, by default 'mag_contour_ccm.pdf'
    """

    
    gonzData = np.loadtxt(datafile,delimiter=',')
    
    # Adjusting the 2MASS Ks magnitude to match the UKIRT K
    # extinction and error values
    gonzData[:,4] = (gonzData[:,4] + 0.689*0.012)/1.069
    gonzData[:,5] = gonzData[:,5]/1.069

    # Adjusting UKIRT H-K color to match 2MASS H-K_s
    # map_data[:,16] = map_data[:,16]*1.062 + 0.017
    # map_data[:,17] = map_data[:,17]*1.062

    if ccm:
        # Calculating the extinction using Cardelli 1989 law
        ak = map_data[:,16] * 1.5
        ak_error = map_data[:,17] * 1.5
    else:
        ak = map_data[:,10] - 13.1
        ak_error = map_data[:,11]



    ukirt_ak = ak
    gonz_ak = gonzData[:,4]
    gonz_error = gonzData[:,5]

    plt.style.use('bmh')
    fig, ax = plt.subplots()
    ax.set_axisbelow(True)
    # ipdb.set_trace()
    ak_df = pd.DataFrame({'UKIRT':ukirt_ak, 'VVV':gonz_ak, 'UKIRT_error':ak_error, 'VVV_error':gonz_error})

    ax.set_axisbelow(True)
    ax.yaxis.grid(color='gray')
    ax.xaxis.grid(color='gray')

    ax.scatter(ukirt_ak,gonz_ak,s=0.1,color='black')
    sns.kdeplot(data=ak_df, x="UKIRT", y="VVV", fill=True, alpha=0.8, levels=10, thresh=0.03, color='blue', ax=ax)
    
    plt.xlabel(r'$A(K)_{UKIRT}$')
    plt.ylabel(r'$A(K)_{VVV}$')

    plt.xlim(0,4)
    plt.ylim(0,4)

    ax.plot([-0.1,4],[-0.1,4], color='red', linestyle='--', linewidth=1.0)

    from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
    from mpl_toolkits.axes_grid1.inset_locator import mark_inset

    x1, x2, y1, y2 = 0.2, .9, 0.2, 0.9 # specify the limits of the zoom
    axins = zoomed_inset_axes(ax, 2.4, loc='upper left') # zoom = 2
    axins.set_axisbelow(True)

    axins.scatter(ukirt_ak,gonz_ak,s=0.1,color='black')
    sns.kdeplot(data=ak_df, x="UKIRT", y="VVV", fill=True, alpha=0.8, levels=10, thresh=0.03, color='blue', ax=axins)
    axins.plot([-0.1,4],[-0.1,4], color='red', linestyle='--', linewidth=1.0)

    axins.set_xlim(x1, x2)
    axins.set_ylim(y1, y2)

    plt.xticks(visible=False)
    plt.yticks(visible=False)
    plt.xlabel('')
    plt.ylabel('')


    plt.tick_params(
        axis='both',
        bottom=False,
        left=False,
        top=False,
        right=False,
        labelbottom=False,
        color='black')

    mark_inset(ax, axins, loc1=3, loc2=4, fc="none", ec="0")

    plt.show()


#------------------------------------------------------------------------------

def plot_ccm_comparison(map_data, path='figs/',figname='mag_contour_ccm.pdf'):
    """Plots the contour map of UKIRT vs VVV extinction values using the Cardelli 1989 law

    Parameters
    ----------
    map_data : np.array
        The map data of the UKIRT data.
    path : str, optional
        The location to save the figure, by default 'figs/'
    figname : str, optional
        The name of the figure including type, by default 'mag_contour_ccm.pdf'
    """
    
    
    # Calculating the extinction using Cardelli 1989 law
    ak_cardelli = map_data[:,16] * 1.5
    ak_error_cardelli = map_data[:,17] * 1.5

    ak_lum_fit = map_data[:,10] - 13.1
    ak_error_lum_fit = map_data[:,11]

    plt.style.use('bmh')
    fig, ax = plt.subplots()
    ax.set_axisbelow(True)
    # ipdb.set_trace()
    ak_df = pd.DataFrame({'Cardelli':ak_cardelli, 'Luminosity':ak_lum_fit, \
            'Cardelli_error':ak_error_cardelli, 'Luminosity_error':ak_error_lum_fit})

    ax.set_axisbelow(True)
    ax.yaxis.grid(color='gray')
    ax.xaxis.grid(color='gray')

    ax.scatter(ak_cardelli,ak_lum_fit,s=0.1,color='black')
    # sns.kdeplot(data=ak_df, x="UKIRT", y="VVV", fill=True, alpha=0.8, levels=10, thresh=0.03, color='blue', ax=ax)
    
    plt.xlabel(r'$A(K)_{Cardelli}$')
    plt.ylabel(r'$A(K)_{Luminosity Fit}$')

    plt.xlim(0,4)
    plt.ylim(0,4)

    ax.plot([-0.1,4],[-0.1,4], color='red', linestyle='--', linewidth=1.0)
    
#------------------------------------------------------------------------------
# Commented below is tools for creating inset plots
#-------------------------------------------------------------------------------

# 
    # from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
    # from mpl_toolkits.axes_grid1.inset_locator import mark_inset
# 
    # x1, x2, y1, y2 = 0.2, .9, 0.2, 0.9 # specify the limits of the zoom
    # axins = zoomed_inset_axes(ax, 2.4, loc='upper left') # zoom = 2
    # axins.set_axisbelow(True)
# 
    # axins.scatter(ukirt_ak,gonz_ak,s=0.1,color='black')
    # sns.kdeplot(data=ak_df, x="UKIRT", y="VVV", fill=True, alpha=0.8, levels=10, thresh=0.03, color='blue', ax=axins)
    # axins.plot([-0.1,4],[-0.1,4], color='red', linestyle='--', linewidth=1.0)
# 
    # axins.set_xlim(x1, x2)
    # axins.set_ylim(y1, y2)
# 
    # plt.xticks(visible=False)
    # plt.yticks(visible=False)
    # plt.xlabel('')
    # plt.ylabel('')
# 
# 
    # plt.tick_params(
        # axis='both',
        # bottom=False,
        # left=False,
        # top=False,
        # right=False,
        # labelbottom=False,
        # color='black')
# 
    # mark_inset(ax, axins, loc1=3, loc2=4, fc="none", ec="0")

    plt.show()
    ipdb.set_trace()
    print(np.max(ak_cardelli-ak_lum_fit))


#------------------------------------------------------------------------------

def reject_outliers(data, m=4):
    return data[abs(data - np.mean(data)) < m * np.std(data)]

#------------------------------------------------------------------------------

if __name__=='__main__':

    # test_map = read_map('maps/map_PSF_2017_2_gonzGrid.map')
    test_map = read_map('maps/map_PSF_2017_2.map')

    # plot_grid_map(test_map,func=rel_diff,axis=10,figname='UKIRTgonzDIFF_CAPPED.pdf',useangle=False)
    import ipdb

    # Histogram differences using normal luminosity fitted extinction
    # df = plot_differences_hist(test_map, func=rel_diff, axis=10)

    # Histogram differences using Cardelli extinction
    # df = plot_differences_hist(test_map, func=rel_diff, axis=16, ccm=True)

    # plot_mags(test_map)
    # plot_mag_card(test_map, cap=True, cap_value=1.5)

    # Creating a contour plot
    # plot_contour(test_map, datafile='../data/gonz/gonzData_2018.dat', ccm=False)
    plot_ccm_comparison(test_map)