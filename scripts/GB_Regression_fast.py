"""
    Perform Gradient Boosting Regression and create a prediction interval
    Used to create a "mean-state: sounding from U of Wyoming Sounding data 
    .csv for this in made in the get_sounding_uw file in the utils folder

    Prediction interval is only used for visual and quality check purposes
    The mean of the prediction interval used will be used to create the wind speed 
    profile

    Following scikit learn's tutorial with a new dataset (nothing fancy done otherwise)
    Link: https://scikit-learn.org/stable/auto_examples/ensemble/plot_gradient_boosting_quantile.html 
    
    Same as GB_Regression but on the top 5% of near surface winds speed days

    October 1, 2023
    lbuchart@eoas.ubc.ca
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import metpy.calc as mpcalc

from metpy.units import units
from sklearn.model_selection import train_test_split
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.metrics import mean_pinball_loss, mean_squared_error
from icecream import ic

from context import utils_dir

##########

my_file = "all_soundings.csv"
df = pd.read_csv(utils_dir + my_file, sep=',')

target = "speed"

# create vector of input values (in this case evenly spaced heights)
def round_down(num, divisor):
    return num - (num%divisor)

# Drop any rows with all NaN values for T, Td, winds
df = df.dropna(
    subset=("temperature", "dewpoint", "direction", "speed", "height"), how="any"
).reset_index(drop=True)

df.drop(list(df.filter(regex = "Unnamed")), axis=1, inplace=True)
ic(df.shape[0])

# assign units to key variables
p = df["pressure"].values * units.hPa  
T = df["temperature"].values * units.degC
Td = df["dewpoint"].values * units.degC 
wind_speed = df["speed"].values * units.knots
wind_dir = df["direction"].values * units.degrees
u, v = mpcalc.wind_components(wind_speed, wind_dir) 
# add the wind components to the dataframe
df["U"] = u
df["V"] = v

print("Get Mixing Ratio")
# get vapour mixing ratio and add to dataframe
rh = mpcalc.relative_humidity_from_dewpoint(T, Td)
Q = mpcalc.mixing_ratio_from_relative_humidity(p, T, rh)
#Q = mpcalc.saturation_mixing_ratio(p, T)  
Q = Q.to("g/kg")
df["mixing_ratio"] = Q

print("Get Potential Temperature")
T = T.to("K")
Th = mpcalc.potential_temperature(p, T)
Th = Th.to("K")
df["PotTemp"] = Th

# instead of elvation ASL, get elevation AGL
df["height"] = df["height"] - df["elevation"][1]


## start modifying the dataframe
# now only get pressures greater than 900 hPa
df_low = df[df["pressure"] > 930]
ic(df_low.shape[0])

# now only get those values which are in the 95th or greater percentile for low level winds
quant = df_low.quantile(q=0.95) 
ic(quant["speed"])
df_low_fast = df_low[df_low["speed"] > quant["speed"]]

# make a final data frame consisting of full sounding data only the high low level winds days
dates = df_low_fast["time"].tolist()
list_dates = [i for n, i in enumerate(dates) if i not in dates[:n]]

mask = df["time"].isin(list_dates)
df_fast = df.loc[mask]

ic(df_fast.shape[0], len(list_dates))

# reorder the data frame according to ascending height (think one very dense sounding)
df_fast = df_fast.sort_values(by="height", ascending=True)

# plot the total wind profile for a sanity check
plt.figure(figsize=(12,12))
plt.plot(df_fast["U"], df_fast["height"], "k.")
plt.title("U Wind Profile")
plt.savefig("fast_u_wind_profile.png")

# plot the total wind profile for a sanity check
plt.figure(figsize=(12,12))
plt.plot(df_fast["V"], df_fast["height"], "k.")
plt.title("V Wind Profile")
plt.savefig("fast_v_wind_profile.png")

## start the gradient boosting regression

# pull out the wind speed variable we are trying to predict
y = np.atleast_2d(df_fast["U"].ravel()).T
X = np.atleast_2d(df_fast["height"]).T

# train test split
random_state = 11
test_size = 0.2
X_train, X_test, y_train, y_test = train_test_split(X, y ,test_size=test_size, random_state=random_state)

all_models = {}
common_params = dict(
    learning_rate = 0.1, 
    n_estimators = 100, 
    max_depth = 3,
    min_samples_leaf = 10, 
    min_samples_split = 10,
)

int_bounds = [0.05, 0.5, 0.90, 0.95]
for alpha in int_bounds:
    gbr = GradientBoostingRegressor(loss="quantile", alpha=alpha, **common_params)
    all_models["q %1.2f" % alpha] = gbr.fit(X_train, y_train)

# also fit a baseline trained with the mean square error
grb_ls = GradientBoostingRegressor(loss="absolute_error", **common_params)
all_models["friedman_mse"] = grb_ls.fit(X_train, y_train)

nhh = round_down( max(df_fast["height"]), 10)
hh = np.arange(0, nhh, 10)
hh = np.append(hh, max(df_fast["height"]))
hh = np.atleast_2d(hh).T

# predict the measured values, 5th, median, and 95th percentiles
y_pred = all_models["friedman_mse"].predict(hh)
y_lower = all_models["q 0.05"].predict(hh)
y_upper = all_models["q 0.95"].predict(hh)
y_med = all_models["q 0.50"].predict(hh)

y_speed = all_models["q 0.90"].predict(hh)

# plot
# plot just the height and the median in the lowest 10000m
fig= plt.figure(figsize=(12, 12))
plt.plot(y_med[0:101], hh[0:101], "r-")
plt.savefig("fast_Low_level_check.png")

# plot everything
fig = plt.figure(figsize=(12, 12))
plt.plot(y_test, X_test, "k.", markersize=10, label="Measured Wind Speed")
plt.plot(y_med, hh, "r-", label="Predicted Median")
plt.plot(y_pred, hh, "b--", label="Predicted Mean")
#plt.plot(y_speed, hh, "y-", linewidth=2, label="LES Sounding")
plt.plot(y_upper, hh, "g-")
plt.plot(y_lower, hh, "g-")
plt.fill_betweenx(hh.ravel(), 
                  y_lower, y_upper, alpha=0.4, facecolor="g",
                  label="Predicted 90% Interval")

plt.title("Observed and Predicted 95th Percentile U Wind Speed - Vernon Sounding")
plt.xlabel("Speed [kts]")
plt.ylabel("Height [m]")
plt.legend()
plt.ylim([0, 10000])
plt.savefig("fast_GBR_Speed.png")


## lets do the same but for the mixing ratio

# plot the total mixing ratio profile for a sanity check
plt.figure(figsize=(12,12))
plt.plot(df_fast["mixing_ratio"], df_fast["height"], "o")
plt.title("Mixing Ratio Profile")
plt.savefig("mixing_ratio_profile.png")

# pull out the mixing ratio variable we are trying to predict
y = np.atleast_2d(df_fast["mixing_ratio"].ravel()).T
X = np.atleast_2d(df_fast["height"]).T

X_train, X_test, y_train, y_test = train_test_split(X, y ,test_size=test_size, random_state=random_state)

all_models = {}
common_params = dict(
    learning_rate = 0.05, 
    n_estimators = 100, 
    max_depth = 2,
    min_samples_leaf = 50, 
    min_samples_split = 50,
)

for alpha in int_bounds:
    gbr = GradientBoostingRegressor(loss="quantile", alpha=alpha, **common_params)
    all_models["q %1.2f" % alpha] = gbr.fit(X_train, y_train)
    
# also fit a baseline trained with the mean square error
grb_ls = GradientBoostingRegressor(loss="squared_error", **common_params)
all_models["mse"] = grb_ls.fit(X_train, y_train)

# plot this thing
# including the measured values, 5th, median, and 95th percentiles
y_pred = all_models["mse"].predict(hh)
y_lower = all_models["q 0.05"].predict(hh)
y_upper = all_models["q 0.95"].predict(hh)
y_med = all_models["q 0.50"].predict(hh)

y_mixing = all_models["q 0.90"].predict(hh)

fig = plt.figure(figsize=(12, 12))
plt.plot(y_test, X_test, "k.", markersize=8, label="Measured Mixing Ratio")
plt.plot(y_med, hh, "r-", label="Predicted Median")
plt.plot(y_pred, hh, "b--", label="Predicted Mean")
plt.plot(y_upper, hh, "g-")
plt.plot(y_lower, hh, "g-")
plt.fill_betweenx(hh.ravel(), 
                  y_lower, y_upper, alpha=0.4, facecolor="g",
                  label="Predicted 90% Interval")

plt.title("Observed and Predicted 95th Percentile Winds Mixing Ratio - Vernon Sounding")
plt.xlabel("Mixing Ratio [g/kg]")
plt.ylabel("Height [m]")
plt.legend()
plt.savefig("fast_GBR_Mixing_Ratio.png")


## lastly do this again with potential temperature 
# (for special experiment only)

# pull out potential temperature variable we are trying to predict
y = np.atleast_2d(df_fast["PotTemp"].ravel()).T
X = np.atleast_2d(df_fast["height"]).T

X_train, X_test, y_train, y_test = train_test_split(X, y ,test_size=test_size, random_state=random_state)

all_models = {}
common_params = dict(
    learning_rate = 0.1, 
    n_estimators = 200, 
    max_depth = 4,
    min_samples_leaf = 50, 
    min_samples_split = 50,
)

for alpha in int_bounds:
    gbr = GradientBoostingRegressor(loss="quantile", alpha=alpha, **common_params)
    all_models["q %1.2f" % alpha] = gbr.fit(X_train, y_train)

# also fit a baseline trained with the mean square error
grb_ls = GradientBoostingRegressor(loss="squared_error", **common_params)
all_models["mse"] = grb_ls.fit(X_train, y_train)

# including the measured values, 5th, median, and 95th percentiles
y_pred = all_models["mse"].predict(hh)
y_lower = all_models["q 0.05"].predict(hh)
y_upper = all_models["q 0.95"].predict(hh)
y_med = all_models["q 0.50"].predict(hh)

y_pottemp = all_models["q 0.90"].predict(hh)

# plot this thing
# including the measured values, 5th, median, and 95th percentile
fig = plt.figure(figsize=(12, 12))
plt.plot(y_test, X_test, "k.", markersize=8, label="Measured Mixing Ratio")
plt.plot(y_med, hh, "r-", label="Predicted Median")
plt.plot(y_pred, hh, "b--", label="Predicted Mean")
plt.plot(y_upper, hh, "g-")
plt.plot(y_lower, hh, "g-")
plt.fill_betweenx(hh.ravel(), 
                  y_lower, y_upper, alpha=0.4, facecolor="g",
                  label="Predicted 80% Interval")

plt.title("Observed and Predicted 95th Percentile Winds Pot Temp - Vernon Sounding")
plt.xlabel("Potential Temperature [K]]")
plt.ylabel("Height [m]")
plt.legend()
plt.savefig("fast_GBR_Pot_Temp.png")

print("Complete")