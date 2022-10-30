#------------------------------------------------------------------------------
# Filename: interpolateMap.py
#
# Given a "raw" map file, interpolates the values of missing or bad pixels.
#
# Programmer: Aiden Zelakiewicz (zelakiewicz.1@osu.edu)
#
# Revision History:
#   06-Oct-2022: File Created
#------------------------------------------------------------------------------

import numpy as np
import scipy.interpolate as interp
import pandas as pd
import ipdb

filename = 'maps/map_PSF_2017_2.map'

map_data = np.loadtxt(filename, delimiter = ',')

# Reduce the map to only parameters we care about
ind = [0,1,10,11,16,17]

extinction_map = map_data[:,ind]

# Find the pixels with no errors
no_error_ak = np.isnan(extinction_map[:,3])
no_error_color = np.isnan(extinction_map[:,5])

# Find indices of pixels with errors within 95%
bad_error_ak = extinction_map[:,3] > np.nanpercentile(extinction_map[:,3], 99.7)
bad_error_color = extinction_map[:,5] > np.nanpercentile(extinction_map[:,5], 99.7)

bad_ind_ak = np.logical_or(no_error_ak, bad_error_ak)
bad_ind_color = np.logical_or(no_error_color, bad_error_color)

bad_map_points_ak = extinction_map[bad_ind_ak][:,0:2]
bad_map_points_color = extinction_map[bad_ind_color][:,0:2]

good_map_values_ak = extinction_map[~bad_ind_ak]
good_map_values_color = extinction_map[~bad_ind_color]
ipdb.set_trace()

# Interpolate the AK and color values
ak_interp = interp.interp2d(good_map_values_ak[:,0], good_map_values_ak[:,1], good_map_values_ak[:,2], bounds_error=True)
color_interp = interp.interp2d(good_map_values_color[:,0], good_map_values_color[:,1], good_map_values_color[:,4], bounds_error=True)

ak_new = np.diagonal(ak_interp(bad_map_points_ak[:,0], bad_map_points_ak[:,1]))
color_new = np.diagonal(color_interp(bad_map_points_color[:,0], bad_map_points_color[:,1]))

new_map_data = map_data
# point_ind_ak = np.logical_and(np.logical_or.reduce([map_data[:,0]==point for point in bad_map_points_ak[:,0]]),
#                                 np.logical_or.reduce([map_data[:,1]==point for point in bad_map_points_ak[:,1]]))

point_ind_ak = np.zeros(len(map_data), dtype=bool)
point_ind_color = np.zeros(len(map_data), dtype=bool)

for (i, pixel) in enumerate(map_data):
    for bad_ak in bad_map_points_ak:
        if np.all(np.equal(pixel[0:2],bad_ak)):
            point_ind_ak[i] = True
            break
    
    for bad_color in bad_map_points_color:
        if np.all(np.equal(pixel[0:2],bad_color)):
            point_ind_color[i] = True
            break

# point_ind_color = np.logical_and([map_data[:,0]==point for point in bad_map_points_color[:,0]],[map_data[:,1]==point for point in bad_map_points_color[:,1]])
# point_ind_color = np.where((map_data[:,0]==bad_map_points_color[:,0]) & (map_data[:,1]==bad_map_points_color[:,1]))
ipdb.set_trace()
new_map_data[:,10][point_ind_ak] = ak_new
new_map_data[:,16][point_ind_color] = color_new

np.savetxt('maps/map_PSF_2017_2_interp.map', new_map_data, delimiter = ',')