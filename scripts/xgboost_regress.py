"""    
    Use extreme gradient boosting to predict a wind speed 
    and mixing ratio sounding
    
    Following tutorial from:
    https://www.datatechnotes.com/2019/06/regression-example-with-xgbregressor-in.html
         
    lbuchart@eoas.ubc.ca
    July 22, 2022
"""

import numpy as np
import pandas as pd
import xgboost as xgb
import matplotlib.pyplot as plt
import metpy.calc as mpcalc

from metpy.units import units
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score, KFold
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import GridSearchCV, RandomizedSearchCV

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

print("Get Mixing Ratio")
# get vapour mixing ratio and add to dataframe
Q = mpcalc.saturation_mixing_ratio(p, T)
Q = Q.to("g/kg")
df["mixing_ratio"] = Q

# instead of elvation ASL, get elevation AGL
df["height"] = df["height"] - df["elevation"][1]
#print(df.describe())
#print("Skew: \n{}".format(df.skew()))

# reorder the data frame according to ascending height (think one very dense sounding)
df = df.sort_values(by="height", ascending=True)
print(df)

# start the Regression
print("Start Regression")
target = "u_wind"

# sample plots
plt.hist(df["mixing_ratio"])
plt.savefig("mixing_hist.png")
plt.close()

plt.hist(df["u_wind"])
plt.savefig("Uwind_hist.png")
plt.close()

outlier_column = "mixing_ratio"
#print( np.log1p(df[outlier_column]).skew() )

# we have shown that the log transform adequately reduces the skewness, now apply it
df[outlier_column] = np.log1p(df[outlier_column])

# ready to start this thing
data_sel = df.copy()

# drop the target and any of the string columns or unnessecary ones (lat, lon, elevation) 
X = data_sel.drop([target, "station", "time", "station_number", "latitude", "longitude", "elevation", "pw"], axis=1)
y = data_sel[target]

# train test split
X_train, X_test, y_train, y_test = train_test_split(X, y ,test_size=0.2, random_state=7, shuffle=True)

# look at decision trees
from sklearn.tree import DecisionTreeRegressor
from sklearn.tree import plot_tree

model_tree = DecisionTreeRegressor(max_depth=8, random_state=44)
model_tree.fit(X_train, y_train)
pred_tree = model_tree.predict(X_test)

# plot the tree
plt.figure(figsize=(20,10))  # set plot size (denoted in inches)
plot_tree(model_tree, fontsize=10, max_depth = 3,feature_names = X.columns, filled = False)
plt.savefig("Decision_Tree.png")

## now actually do the XGBoost Stuff ##
# train test split
random_state = 11
test_size = 0.2
X_train, X_test, y_train, y_test = train_test_split(X, y ,test_size=test_size, random_state=random_state, shuffle=True)

n_estimators = 100
max_depth = 5
reg_alpha = 0.5
#model_xgb = xgb.XGBRegressor(n_estimators=n_estimators, 
#                             max_depth=max_depth, 
#                             reg_alpha=reg_alpha, 
#                             random_state=random_state)
 
model_xgb = xgb.XGBRegressor(base_score=0.5, booster='gbtree', colsample_bylevel=1,
       colsample_bynode=1, colsample_bytree=1, gamma=0,
       importance_type='gain', learning_rate=0.1, max_delta_step=0,
       max_depth=2, min_child_weight=1, missing=None, n_estimators=100,
       n_jobs=1, nthread=None, objective='reg:linear', random_state=random_state,
       reg_alpha=20, reg_lambda=1, scale_pos_weight=1, seed=None,
       silent=None, subsample=1, verbosity=1)

eval_set = [(X_train, y_train), (X_test, y_test)]
model_xgb.fit(X_train, y_train, eval_metric=["rmse"], eval_set=eval_set, verbose=True)
preds_xgb1 = model_xgb.predict(X_test)

rmse = np.sqrt(mean_squared_error(np.expm1(y_test), np.expm1(preds_xgb1)))
#print("The new RMSE: %f" % (rmse))
#print( model_xgb.get_xgb_params() )

#print("Performance metrics")
# Retrieve performance metrics
results = model_xgb.evals_result()
epochs = len(results['validation_0']['rmse'])
x_axis = range(0, epochs)
# plot RMSE
fig, ax = plt.subplots()
ax.plot(x_axis, results['validation_0']['rmse'], label='Train')
ax.plot(x_axis, results['validation_1']['rmse'], label='Test')
ax.legend()
plt.ylabel('RMSE')
plt.title('XGBoost RMSE')
plt.savefig("XGBoost_performance.png")

# plot the values
plt.figure(figsize=(12, 12))
x_ax = range(len(y_test))
plt.plot(y_test, X_test["height"], "-o", label="Original")
plt.plot(preds_xgb1, X_test["height"], "-o", label="Predicted")
plt.legend()
plt.title("U_wind Comparison")
plt.savefig("U_wind_comparison.png")

# now use the XGBoost Grid Search Feature
xgb1 = xgb.XGBRegressor() 

parameters = {
            'objective': ['reg:squarederror'],
            'booster':['gbtree'], # ,'gblinear'
            'max_depth': [2,4,6],
            'learning_rate': [0.1,0.5, 1],
            'n_estimators': [100, 200],
            'colsample_bytree': [0.3, 0.7],
            'importance_type': ['gain'],
            'subsample': [0.5, 1],
            'reg_alpha': [0.5,1, 5],
            'reg_lambda': [2,5]}

xgb_grid = GridSearchCV(xgb1,
                        parameters,
                        cv = 2,
                        verbose=2,
                        scoring='neg_root_mean_squared_error') # https://scikit-learn.org/stable/modules/model_evaluation.html

xgb_grid.fit(X_train,y_train)

print(xgb_grid.best_score_)
print(xgb_grid.best_params_)
print(xgb_grid.best_estimator_)

xg_reg = xgb_grid.best_estimator_
eval_set = [(X_train, y_train), (X_test, y_test)]
xg_reg.fit(X_train, y_train, eval_metric=["rmse"], eval_set=eval_set, verbose=True)
preds = xg_reg.predict(X_test)

rmse = np.sqrt(mean_squared_error(np.expm1(y_test), np.expm1(preds)))
print("RMSE: %f" % (rmse))

esults = model_xgb.evals_result()
epochs = len(results['validation_0']['rmse'])
x_axis = range(0, epochs)
# plot RMSE
fig, ax = plt.subplots()
ax.plot(x_axis, results['validation_0']['rmse'], label='Train')
ax.plot(x_axis, results['validation_1']['rmse'], label='Test')
ax.legend()
plt.ylabel('RMSE')
plt.title('XGBoost RMSE')
plt.savefig("SEARCH_XGBoost_performance.png")

# plot the values
plt.figure(figsize=(12, 12))
x_ax = range(len(y_test))
plt.plot(y_test, X_test["height"], "-o", label="Original")
plt.plot(preds_xgb1, X_test["height"], "-o", label="Predicted")
plt.legend()
plt.title("U_wind Comparison")
plt.savefig("SEARCH_U_wind_comparison.png")

plt.figure(figsize=(12, 12))
plt.plot(y_train, X_train["height"], "-o", label="Training")
plt.plot(preds_xgb1, X_test["height"], "-o", label="Predicted")
plt.legend()
plt.title("U Winds Train")
plt.savefig("Training_U_Wind")