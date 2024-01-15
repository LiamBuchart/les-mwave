import context
import numpy as np
import pandas as pd
from datetime import datetime, timedelta
import matplotlib.pyplot as plt

import metpy.calc as mpcalc
from metpy.units import units
from siphon.simplewebservice.wyoming import WyomingUpperAir

from context import root_dir, utils_dir

##########
"""
    Get upper air data for a specified station from the
    University of Wyoming Sounding Page
    https://weather.uwyo.edu/upperair/sounding.html
    
    Usage: Generate a mean summer wind profile to be used 
    for a WRF-LES simulation
    
    Modified script - original author: Chris Rodell
    https://github.com/cerodell
    
    contact: lbuchart@eoas.ubc.ca
"""
pd.get_option('display.max_rows', None)
    
station = "CWVK"
s_elev = 379  # station elevation at Vernon
start_date = datetime(2018, 5, 1, 12)  # start of 2018 fire season
end_date = datetime.today()
time_step = 1  # iterate every day
time = 00  # which sounding do we want: 12Z

mstart = 5  # defining the month in which the fire season starts
mend = 9  # defining the month in which the fire season ends

delta = end_date - start_date
print("Number of Days: ", delta.days)

###########

def daterange(dstart, dend):
    # create a list of the dates between start and end spacing one day
    delta = dend - dstart
    diff = delta.days
    for n in range(int(diff)):
        yield dstart + timedelta(n)

###########

# check and see if we have downloaded any soundings before if not create a dataframe for them
# or can be used similarly if you move the all soundings file out and then can create new
try:   
    all_soundings = pd.read_csv(utils_dir + "all_soundings.csv", sep=",")
    print("file exists - openings")
except IOError:
    print("file does not exist or is not accesible, create")
    all_soundings = pd.DataFrame()  # empty datframe to store all data in

print(all_soundings)

result = all_soundings.dtypes   #isin([str(start_date)]).any().any()
print("Output: ")
print(type(str(start_date)), " ", result)

for date in daterange(start_date, end_date):
    
    # loop through dates and append sounding info to the all_soundings dataframess
    month = str(date)[5:7]
    yr = str(date)[0:4]
    day = str(date)[8:10]
    duse = int( day+month+yr )
    month = int(month)
    
    #print(duse, " is it in here? ", duse in all_soundings.values)

    if month < mstart or month > mend or duse in all_soundings.values:
        # dont want to make our wind profile from outside of fire season
        print("Already there or out of season")
        pass
    else:  # get months inside the fire season
        try: 
            df = WyomingUpperAir.request_data(date, station)

            df["ddmmyyyy"] = duse
            all_soundings = pd.concat([all_soundings, df])
        except:
            print("No Sounding Data Today")
            pass

        print(str(date), " - ", len(all_soundings["height"]))
        
        all_soundings.to_csv("all_soundings.csv", sep=',')

print("Complete")
