# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 12:19:55 2016

Concatenate files for velocity data

@author: hkarimi
"""
import matplotlib.pyplot as plt
import numpy as np
from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from sklearn.svm import SVR
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.model_selection import GridSearchCV
from sklearn.base import clone

# load the training data, coords
coords = np.loadtxt('./concatenated/coords.txt')
    
# load the targets
ux = np.loadtxt('./concatenated/ux.txt')
uy = np.loadtxt('./concatenated/uy.txt')
uz = np.loadtxt('./concatenated/uz.txt')

## split into test and train data
coords_train, coords_test, ux_train, ux_test = train_test_split(
         coords, ux, test_size=0.1, random_state=0)
coords_train, coords_test, uy_train, uy_test = train_test_split(
         coords, uy, test_size=0.1, random_state=0)
coords_train, coords_test, uz_train, uz_test = train_test_split(
         coords, uz, test_size=0.1, random_state=0)
#%% train by MLPRegressor 
# build a pipe-line to standardize the input coordinates
pipe_ux_mlpr = Pipeline([ ('scl', StandardScaler()),
                     ('regr',MLPRegressor(activation='tanh')) ])
learning_range = [0.0001, 0.001, 0.01, 0.1, 1]
layer_range = [ (10,5), (10,10), (10,10,5)]
param_grid = [{'regr__learning_rate_init': learning_range,
              'regr__hidden_layer_sizes': layer_range}]
                     
gs = GridSearchCV(estimator=pipe_ux_mlpr,
                  param_grid=param_grid,
                  verbose=1,
                  cv=3,
                  n_jobs=2)
                  
gs.fit(coords_train,ux_train)
print(gs.best_score_)
print(gs.best_params_)

#%% use the best parameters from the block above:
pipe_ux_mlpr = Pipeline([ ('scl', StandardScaler()),
                          ('regr',MLPRegressor(learning_rate_init=0.1,
                                          hidden_layer_sizes=(10,10,5),
                                          activation='tanh'))
                                          ])
pipe_uz_mlpr = clone(pipe_ux_mlpr)

pipe_ux_mlpr.fit(coords_train,ux_train)
pipe_uz_mlpr.fit(coords_train,uz_train)
print('Train Accuracy: %.3f' % pipe_ux_mlpr.score(coords_train, ux_train))
print('Test Accuracy: %.3f' % pipe_ux_mlpr.score(coords_test, ux_test))


#%% train by SVR
#pipe_ux_svr = Pipeline([ ('scl',StandardScaler()),
#                         ('clf',SVR(kernel='rbf')) ])
#
#param_range = [0.001, 0.01, 0.1, 1.0, 10.0, 100.0, 1000.0]  
#param_grid = [{'clf__C': param_range, 
#               'clf__gamma': param_range}]
#                  
#gs = GridSearchCV(estimator=pipe_ux_svr,
#                  param_grid=param_grid,
##                  scoring='accuracy',
#                  verbose=1,
#                  cv=3,
#                  n_jobs=1)
                  
#gs.fit(coords_train,ux_train)
#print(gs.best_score_)
#print(gs.best_params_)
# Note: best_params_ = {'clf__gamma': 1.0, 'clf__C': 0.1}


#%% use the best parameters from the block above:
pipe_ux_svr = Pipeline([ ('scl',StandardScaler()),
                         ('clf',SVR(kernel='rbf',
                                    C=0.1,
                                    gamma=1.0))
                                    ])
pipe_ux_svr.fit(coords_train,ux_train)
#print('Train Accuracy: %.3f' % pipe_ux_svr.score(coords_train, ux_train))
#print('Test Accuracy: %.3f' % pipe_ux_svr.score(coords_test, ux_test))


#%%
## test predictions
# build test coordinates
phi = 0.0; z = -0.0; b = 0.8; sd = 2;
testCoords = np.array([]).reshape(0,5)
R = np.linspace(0,1,20)
for r in R:
    coord = [[r, phi, z, b, sd]]
    testCoords = np.append(testCoords, coord, axis = 0)
    
uzPred = pipe_ux_mlpr.predict(testCoords)
#uxPred = pipe_ux_svr.predict(testCoords)
f, ax = plt.subplots()
ax.xaxis.get_major_formatter().set_powerlimits((0,1))
ax.plot(uzPred,R)
plt.show()
#%%

#pipe_svc = Pipeline([('scl', StandardScaler()),
#            ('clf', SVC(random_state=1))])
#
#gs = GridSearchCV(estimator=pipe_svc, 
#                  param_grid=param_grid, 
#                  scoring='accuracy', 
#                  cv=10,
#                  n_jobs=-1)
#gs = gs.fit(X_train, y_train)
#print(gs.best_score_)
#print(gs.best_params_)