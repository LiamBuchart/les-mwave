"""
 
 Function to plot the scorer parameter l^2(z) = N^2/U^2
 Note that the second (d^2U/dz^2)/U is neglected in this case
 Similar to [enter citation on Hong Kong Gravity Waves]
 
 lbuchart@eoas.ubc.ca
 Created: August 17, 2017   
    
"""

import os
import numpy as np
import metpy.calc as mpcalc
import matplotlib.pyplot as plt

from netCDF4 import Dataset
from context import name_dir
from file_funcs import (setup_script, get_heights, dist_from_ridge)
from scipy.interpolate import UnivariateSpline

from wrf import (getvar, extract_times, 
                 xy, interp2dxy, to_np,
                 interpline, CoordPair)

##########

## USER INPUTS ##
# options ["real", "stable", "neutstab", "logT", "high_step", "low_step"]
exp = "high_step/"
# cross section coordinates
start = (0, 100)
end = (-1, 100)
# where we will get the scorer profile cross section spanning [start, end]
profile_point = 200

## END USER INPUTS ##

# get files and paths set up
path, save_path, relevant_files, wrfin = setup_script(exp)

scorer_profile = []
times = []
# loop through all files and create plots
for ii in range(0, len(wrfin)):
    ncfile = wrfin[ii]

    # get the time in datetime format
    ct = extract_times(ncfile, timeidx=0)
    times.append(str(ct)[0:19])
    print(ct)

    # extract the buoyancy frequency and heights
    U = getvar(ncfile, "ua", units="m s-1", timeidx=0, meta=True)
    theta = getvar(ncfile, "theta", units="K")
    heights = getvar(ncfile, "height", meta=True)
    
    N = mpcalc.brunt_vaisala_frequency(heights, theta)

    # take a cross section of buoyancy frequency, U, and heights through middle of domain along length
    U_line = xy(U, start_point=start, end_point=end)
    U_cross = interp2dxy(U, U_line)
    print("Wind is processed")

    N_line = xy(N, start_point=start, end_point=end)
    N_cross = interp2dxy(N, N_line)
    print("Buoyancy is processed")

    H_line = xy(heights, start_point=start, end_point=end)
    H_cross = interp2dxy(heights, H_line)
    print("Heights are processed")

    # calculate the first term of the scorer parameter (N^2/U^2)
    first_term = (np.square(to_np(N_cross))) / (np.square(to_np(U_cross)))

    print("The mean is: ", np.nanmean(first_term), np.mean(first_term))

    # calculate the second term
    # have to take the second derivative of U by taking the spline of (U(z))
    x = H_cross[:, 0]
    y = U_cross[:, 0]

    # take a cubic spline
    y_spl = UnivariateSpline(x, y, s=0, k=3)

    # quick plot of our spline to ensure things look nice
    plt.figure(figsize=(12,12))
    plt.plot(y, x, 'ro', label="Sounding")
    x_range = np.linspace(x[0], x[-1], 1000)
    plt.plot(y_spl(x_range), x_range, label="Cubic Spline")
    plt.legend()
    plt.xlabel("U Wind [m/s]", fontsize=12)
    plt.ylabel("Height [m]", fontsize=12)
    plt.title("Cubic Spline Check")
    plt.savefig(save_path + "spline_check")

    # take the second derivative
    y_spl_2d = y_spl.derivative(n=2)

    # quick plot of our second derivative to ensure things look nice
    plt.figure(figsize=(12,12))
    plt.plot(y_spl_2d(x_range), x_range, 'k-', label="Second Derivative")
    x_range = np.linspace(x[0], x[-1], 1000)
    plt.legend()
    plt.ylabel("Height [m]", fontsize=12)
    plt.title("Second Derivative Check")
    plt.savefig(save_path + "second_derivative_check")

    # second derivative for the z values we have
    d2udz2 = y_spl_2d(x)

    # get the second term of the scorer parameter in same manner put loop entire array
    shape_xs = np.shape(H_cross)
    sec_der = np.zeros((shape_xs))
    for jj in range(0, shape_xs[1]):
        x = H_cross[:, jj]
        y = U_cross[:, jj]

        y_spl = UnivariateSpline(x, y, s=0, k=3)
        y_spl_2d = y_spl.derivative(n=2)
        d2udz2 = y_spl_2d(x)
        sec_der[:, jj] = d2udz2

    second_term = sec_der / U_cross

    print(np.nanmean(second_term))

    # the scorer parameter
    l2 = first_term - second_term

    # terrain
    ter = getvar(ncfile, "ter", meta=True)
    ter_cross = interpline(ter, 
                           start_point=CoordPair(start[0], start[1]),
                           end_point=CoordPair(end[0], end[1]))
    ridge_dist = dist_from_ridge(ter_cross)

    # plot a contour image of the scorer and a vertical profile in the lee of the ridge
    hh = H_cross[:, 0]
    # make the figure
    fig, axs = plt.subplots(2, 1, 
                            sharex=True,
                            gridspec_kw={"height_ratios":[5, 1]})

    scorer_levels = np.arange(-0.0001, 0.0001, 0.000001)  # contour lines 

    scorer_contour = axs[0].contourf(ridge_dist,
                                hh,
                                to_np(l2),
                                cmap="seismic",
                                extend="both", 
                                levels=scorer_levels)

    ht_fill = axs[1].fill_between(ridge_dist, 0, to_np(ter_cross),
                                  facecolor="saddlebrown")
    
    # hide x labels on all but the bottom 
    for ax in axs:
        ax.label_outer()
        
    # beautify
    axs[0].set_yticks(np.arange(0, 11000, 500))
    axs[0].set_ylim([0, 2500])
    
    axs[1].set_yticks(np.arange(0, 1250, 250))
    axs[1].set_ylim([0, 800])

    # colorbar 
    fig.tight_layout()  # call this before calling the colorbar and after calling 
    cbar = plt.colorbar(scorer_contour, ax=axs)
    cbar.set_label("Scorer Parameter [$m^{-2}$]", fontsize=10)

    plt.savefig(save_path + "scorer_parameter_" + str(ct)[11:19] + "_contour_xs")
    plt.close()
    
    # save a profile at each time
    scorer_profile.append(to_np(l2[:, 4]))

# plot an upstream (close to the edge) scorer parameter profile
shape_profile = np.shape(scorer_profile)
plt.figure(figsize=(12,12))
for tt in range(0, shape_profile[0]):
    plt.plot(scorer_profile[tt], H_cross[:, 5], "-o", label=times[tt], linewidth=2)
 
plt.legend()   
plt.xlabel("Scorer Parameter $l^{2}$ $[m^{-2}]$", fontsize=14)
plt.ylabel("Height AGL [m]", fontsize=14)

plt.ylim([0, 2000])

plt.savefig(save_path + "scorer_profile_all_times")

print("Complete")