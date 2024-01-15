"""
    
   Make nice potential temperature plots 
    
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import json

# from GB_Regression import y_mixing, y_speed, hh
from context import script_dir
from GB_Regression_fast import y_mixing, y_speed, y_pottemp, hh
from get_heights import heights
from file_funcs import extract_val, check_folder
from icecream import ic

from context import json_dir

###########
exps =  ["stable", "neutstab", "log", "low_step", "high_step"]
df = pd.DataFrame(columns=exps)

# folder to save figures (check its there, if not make)
fig_folder = check_folder(script_dir, "figures/paperplots/")

# get mountian height from our json dictionary
with open(str(json_dir) + "config.json") as f:
    config = json.load(f)
m_height = config["topo"]["mht"]  # topography height
T0 = config["sounding"]["T_sfc"]  # surface temperature
P0 = config["sounding"]["P_sfc"]  # surface pressure

for theta_type in exps:

    print(theta_type)
    input_sounding = pd.DataFrame(columns=["height", "theta", "q", "U", "V"])

    hh = hh.ravel()
    ic(hh)
    ind = np.where(hh == 0)
    qfind = y_mixing[ind]

    q1 = round(qfind[0])

    # overwrite the base input_sounding dateframe to account for first row formatting
    # note that the first row of height is actually the surface pressure
    weather_data = {"height":[P0], "theta":[T0], "q":[int(q1)], "U":[np.nan], "V":[np.nan]}
    input_sounding = pd.DataFrame(weather_data, index=[0])

    # build potential temperature profile
    theta = np.zeros((len(hh),), dtype=np.float32)
    if theta_type == "neutstab":

        # neutral lower layer then slowly increaing potential temp
        ind = np.where(hh == int(round((1.5*m_height/10))*10)) # round(number/10)*10
        for ii in range(0, int(ind[0])):
            theta[ii] = T0
    
        slope = 0.005  # potential temp increase / m (theta = T0 + slope*(hh[ii]-hh[ind]))
        for ii in range((int(ind[0])), len(hh)):
            h_change = hh[int(ind[0])]
            theta[ii] = round( T0 + ( slope * (hh[ii] - h_change) ), 2 )

    elif theta_type == "stable":

        # stable and slowing increasing potential temperature with height
        slope = 0.005

        for ii in range(0, len(hh)):
            theta[ii] = round( T0 + ( slope * hh[ii] ), 2 )

    elif theta_type == "log":

        # a logarithmic potential temperature profile (theta = amp*log(str*hh[ii]))
        amp = 3
        str = 2
        for ii in range(0, len(hh)):
            #theta[ii] = T0 + np.sqrt(hh[ii]/1000)
            theta[ii] = T0 + (amp * np.log10(str * hh[ii]))
        
    elif theta_type == "low_step" or theta_type == "high_step":
    
        # sharp potential temp inversion below mountaintop height
        if theta_type == "low_step":
            inv_height = m_height * 0.75  # ensure this value is a multiple of 10 (in hh)
        elif theta_type == "high_step":
            inv_height = m_height * 1.25
        inv_step = 8  # size of the potential temp inversion
    
        ind = np.where(hh == (round(inv_height/10)*10))  # round to nearest 10
    
        count = 0
        len_inv = 160
        for ii in range(0, int(ind[0])):
            theta[ii] = T0     
        for ii in range(int(ind[0]), int(ind[0])+len_inv):
            count = count + 1
            theta[ii] = T0 + (inv_step * (count / len_inv))
            #theta[ii] = (T0 + (T0 + inv_step)) / 2  # intermediate step to make inversion less steep and more physical 
        for ii in range(int(ind[0])+len_inv, len(hh)):  
            theta[ii] = T0 + inv_step
        
    elif theta_type == "real":
    
        theta = y_pottemp
    
    else:
    
        print("No sounding found")
        print("Check Option input line 39")
        quit
    
    # trim to just get the lowest 10 000m
    theta = np.where(hh < 10000, theta, -10)
    hh = np.where(hh < 10000, hh, -10)
    
    tt = []
    hhh = []
    real = []
    uu = []
    qq = []
    for ii in range(len(hh)):
        if hh[ii] != -10:
            tt.append(theta[ii])
            hhh.append(hh[ii])
            real.append(y_pottemp[ii])
            uu.append(y_speed[ii])
            qq.append(y_mixing[ii])
    
    print("Quick peak: ", tt)
    df[theta_type] = tt
    
    #### Do the sam for the real experiment!!!!

ic(df.head())    
heights = heights[1:]  # recall that the 0 index in heights is actually the sfc pressure
    
# plt that single theta profile
rows = 3
cols = 2
fig, axs = plt.subplots(rows, cols, sharex=True)

axs = axs.ravel()

for ii in range(6):

    axs[ii].set_yticks((np.arange(0, 12500, 2500)))
    axs[ii].set_yticklabels(np.arange(0, 12500, 2500), rotation=45)
    if ii == 5:
        theta = y_pottemp
        
    if ii == 5:
        axs[ii].plot(real, hhh, color='k')
        axs[ii].text(0.50, 0.40, "real", transform=axs[ii].transAxes, fontsize=14,
                     verticalalignment='top')
    else:
        axs[ii].plot(df[exps[ii]], hhh, color='k')
        axs[ii].text(0.50, 0.40, exps[ii], transform=axs[ii].transAxes, fontsize=14,
                     verticalalignment='top')
 
axs[0].set_ylabel("Height [m]", fontsize=12)   
axs[2].set_ylabel("Height [m]", fontsize=12) 
axs[4].set_ylabel("Height [m]", fontsize=12) 
axs[-2].set_xlabel("Potential Temperature [K]", fontsize=12)
axs[-1].set_xlabel("Potential Temperature [K]", fontsize=12)
plt.savefig(fig_folder + "/all_Theta_Profile.png")

# plot the speed in a similar manner
fig, ax = plt.subplots(1, 1)
ax.plot(uu, hhh, color='k')
ax.set_yticks((np.arange(0, 12500, 2500)))
ax.set_yticklabels(np.arange(0, 12500, 2500), rotation=45)
ax.set_ylabel("Height [m]", fontsize=12)
ax.set_xlabel("Horizontal Speed [m/s]", fontsize=12)
plt.savefig(fig_folder + "/figure1_speed.png")

# plot the mixing ratio in a similar manner
fig, ax = plt.subplots(1, 1)
ax.plot(qq, hhh, color='k')
ax.set_yticks((np.arange(0, 12500, 2500)))
ax.set_yticklabels(np.arange(0, 12500, 2500), rotation=45)
ax.set_ylabel("Height [m]", fontsize=12)
ax.set_xlabel("Mixing Ratio [g/kg]", fontsize=12)
plt.savefig(fig_folder + "/figure1_mixing.png")

