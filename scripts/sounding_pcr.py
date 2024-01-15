"""
    
    Function which uses Principal Component Regression
    https://towardsdatascience.com/principal-component-regression-clearly-explained-and-implemented-608471530a2f
    to create vertical wind and mixing ratio profile to be used to force wrf-les simulations
        
    contact: lbuchart@eoas.ubc.ca
    
"""

import numpy as np
import pandas as pd
import metpy.calc as mpcalc
from metpy.units import units
import matplotlib.pyplot as plt

from sklearn.preprocessing import scale 
from sklearn.linear_model import LinearRegression, RidgeCV, LassoCV
from sklearn.model_selection import KFold, cross_val_score, train_test_split
from sklearn.metrics import mean_squared_error
from sklearn.decomposition import PCA

from context import utils_dir

##########

my_file = "all_soundings.csv"
df = pd.read_csv(utils_dir + my_file, sep=',')

# Drop any rows with all NaN values for T, Td, winds
df = df.dropna(
    subset=("temperature", "dewpoint", "direction", "speed"), how="all"
).reset_index(drop=True)
 
# assign units to key variables
p = df["pressure"].values * units.hPa  
T = df["temperature"].values * units.degC
Td = df["dewpoint"].values * units.degC 
wind_speed = df["speed"].values * units.knots
wind_dir = df["direction"].values * units.degrees
u, v = mpcalc.wind_components(wind_speed, wind_dir) 

print("get mixing ratio")
# get vapour mixing ratio and add to dataframe
Q = mpcalc.saturation_mixing_ratio(p, T)
Q = Q.to('g/kg')
df["mixing_ratio"] = Q

# instead of elvation ASL, get elevation AGL
df["height"] = df["height"] - df["elevation"][1]

print(df.head())
# start pcr analysis ########## 
print("Starting PCR analysis")

target = "u_wind"
#target2 = "v_wind"
#target = "mixing_ratio"

# drop the target and any of the string columns or unnessecary ones (lat, lon, elevation) 
X = df.drop([target, "station", "time", "station_number", "latitude", "longitude", "elevation", "pw"], axis=1)
y = df[target]

# divide in to train and test values
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# standardize/scale
X_train_scaled = scale(X_train)
X_test_scaled = scale(X_test)

# run baseline regression model (standard linear regression)
cv = KFold(n_splits=10, shuffle=True, random_state=42)

lin_reg = LinearRegression().fit(X_train_scaled, y_train)
lr_score_train = -1 * cross_val_score(lin_reg, X_train_scaled, y_train, cv=cv, scoring='neg_root_mean_squared_error').mean()
lr_score_test = mean_squared_error(y_test, lin_reg.predict(X_test_scaled), squared=False)

print(lr_score_test)

# generate principal components
pca = PCA()  
X_train_pc = pca.fit_transform(X_train_scaled)

# check the principal components
#print( pd.DataFrame(pca.components_.T).loc[:4, :] )

# view the explained variance ratio
print(pca.explained_variance_ratio_)

# determine the number of principal components
lin_reg = LinearRegression()

# empty list to store rmse
rmse_list = []

# loop through principal components
for ii in range(1, X_train_pc.shape[1]+1):
    rmse_score = -1 * cross_val_score(lin_reg, 
                                      X_train_pc[:, :ii], 
                                      y_train, 
                                      cv=cv, 
                                      scoring="neg_root_mean_squared_error").mean()
    
    rmse_list.append(rmse_score)
    
# plot RMSE vs number of principal components
plt.plot(rmse_list, '-o')
plt.xlabel("Number of Principal Components in Regression")
plt.ylabel("RMSE")
plt.title("u winds")
plt.savefig("rmse_pca.png")
plt.close()

# run pcr with chosen principal components 
best_pc_num = 6

# train model 
lin_reg_pc = LinearRegression().fit(X_train_pc[:, :best_pc_num], y_train)

# get cross-validation RMSE (train set)
pcr_score_train = -1 * cross_val_score(lin_reg_pc, 
                                       X_train_pc[:, :best_pc_num],
                                       y_train,
                                       cv=cv,
                                       scoring="neg_root_mean_squared_error").mean()

# train model on training set
lin_reg_pc = LinearRegression().fit(X_train_pc[:, :best_pc_num], y_train)

# get the first principal components that you want
X_test_pc = pca.transform(X_test_scaled)[:, :best_pc_num]

# predict of test data
preds = lin_reg_pc.predict(X_test_pc)

print(len(preds), len(df["height"]))

# plt the predictions
plt.plot(df["u_wind"], df["height"], "o", color="blue")
#plt.plot(preds, "-o", color="red")
plt.xlabel("wind speed")
plt.ylabel("height")
plt.title("Wind Profile")
plt.xlim((-30, 30))
plt.savefig("wind_profile.png")
plt.close()

plt.plot(preds, "-o", color="red")
plt.xlim((0, len(preds)))
plt.savefig("test_set.png")