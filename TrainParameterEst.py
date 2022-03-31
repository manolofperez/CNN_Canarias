import sys, os
import numpy as np
import time
import random
from random import shuffle, choice
import keras
import tensorflow as tf
from tensorflow.keras.models import Model
from tensorflow.keras.regularizers import l2
from tensorflow.keras.constraints import max_norm
from tensorflow.keras.callbacks import EarlyStopping
from tensorflow.keras import backend as K
from tensorflow.keras.layers import *
from tensorflow.keras.optimizers import *
from tensorflow.keras.models import load_model
from keras.utils import np_utils

batch_size = 500
epochs = 500

def readDemogParams(demogParamPath):
    params = []
    first = True
    with open(demogParamPath) as demogParamFile:
        for line in demogParamFile:
            params.append([float(x) for x in line.strip().split()])
    return params
    
def create_cnn(xtest):
	inputShape = (imgRows, imgCols)
	inputs = Input(shape=inputShape)
	x = inputs
	x = Conv1D(125, kernel_size=3, activation='relu', use_bias=False, input_shape=(imgRows, imgCols))(x)
	x = BatchNormalization()(x)
	x = Conv1D(250, kernel_size=3, use_bias=False, activation='relu')(x)
	x = BatchNormalization()(x)
	x = Conv1D(250, kernel_size=3, use_bias=False, activation='relu')(x)
	x = BatchNormalization()(x)
	x = MaxPooling1D(pool_size=3)(x)
	x = Flatten()(x)
	x = Dense(125, activation='relu')(x)
	x = Dropout(0.5)(x)
	x = Dense(125, activation='relu')(x)
	x = Dropout(0.5)(x)
  	# The final fully-connected layer will have a dense layer
	x = Dense(numParams)(x)
	# Construct the CNN
	model = Model(inputs, x)
	# Return the CNN
	return model

x = np.load("../Mod_SSHBackCol.npy",mmap_mode='r')


x=np.array(x[:,0:1000,:])
#transform major alleles in -1 and minor 1
for arr,array in enumerate(x):
  for idx,row in enumerate(array):
    if np.count_nonzero(row) > len(row)/2:
      x[arr][idx][x[arr][idx] == 1] = -1
      x[arr][idx][x[arr][idx] == 0] = 1
    else:
      x[arr][idx][x[arr][idx] == 0] = -1

imgRows, imgCols = x.shape[1:]

demogParams = readDemogParams('../parSSHBackCol.txt')

y = np.array(demogParams)
numParams=y.shape[1]

del (demogParams)

shf = list(range(len(x)))
shuffle(shf)

y = y[shf]
x = x[shf]

#Add missing data as 0s, according to a specifies missing data percentage
#584,872 SNPs and 25,867 missing genotypes = 4.4%
missD_perc = 4.4
missD = int(x.shape[1]*x.shape[2]*(missD_perc/100))
for i in range(x.shape[0]):
  for m in range(missD):
    j = random.randint(0, x.shape[1] - 1)
    k = random.randint(0, x.shape[2] - 1)
    x[i][j][k] = 0

del(missD)

yMeans=np.mean(y, axis=0)
yStds=np.std(y, axis=0)
y = (y-yMeans)/yStds

print (yMeans)
print (yStds)

# Separate train (75%) and validate (25%) sets.
xtrain, xtest = x[int(len(y)*.25):], x[:int(len(y)*.25)]
ytrain, ytest = y[int(len(y)*.25):], y[:int(len(y)*.25)]
del(x)

# Create the CNN network
cnn = create_cnn(xtest)

opt = SGD(lr=0.001)
#opt = 'Adam'
cnn.compile(loss='mean_squared_error', optimizer=opt)

print(cnn.summary())

earlyStopping = EarlyStopping(monitor='val_loss', patience=150, verbose=0, mode='min', restore_best_weights=True)

start = time.time()
cnn.fit(xtrain, ytrain, batch_size=batch_size,
          epochs=epochs,
          verbose=1,
          validation_data=(xtest, ytest),callbacks=[earlyStopping])
print (f'Time: {time.time() - start}')

with open('Trained_Params_10KSims.acc.mod', "w+") as modFile:
    modFile.write(cnn.to_json())
    
# Load the simulations.
x_test = np.load("../testSims/Mod_SSHBackCol.npy",mmap_mode='r')
x_test = np.array(x_test[:,0:1000,:])

# Convert the reference allele to -1.
#transform major alleles in -1 and minor 1
for arr,array in enumerate(x_test):
  for idx,row in enumerate(array):
    if np.count_nonzero(row) > len(row)/2:
      x_test[arr][idx][x_test[arr][idx] == 1] = -1
      x_test[arr][idx][x_test[arr][idx] == 0] = 1
    else:
      x_test[arr][idx][x_test[arr][idx] == 0] = -1

#Add missing data as 0s, according to a specifies missing data percentage
#584,872 SNPs and 25,867 missing genotypes = 4.4%
missD_perc = 4.4
missD = int(x_test.shape[1]*x_test.shape[2]*(missD_perc/100))
for i in range(x_test.shape[0]):
  for m in range(missD):
    j = random.randint(0, x_test.shape[1] - 1)
    k = random.randint(0, x_test.shape[2] - 1)
    x_test[i][j][k] = 0

# Predict parameters for each simulation.
pred = cnn.predict(x_test)

# Save the obtained predictions.
np.savetxt("testSet_ParameterPredictions.txt", pred)

# Load empirical data.
infile=np.loadtxt("../input_Euphorbia.txt")
inp=np.array(infile)

# Create 100 subsets containing 1,000 random SNPs from the full empirical data.
num_samples=100
res = []
for i in range(0,num_samples):
	idx = np.random.choice(inp.shape[0], 1000, replace=False)
	n = inp[idx,:]
	res.append(np.array(n))

# Predict parameters.
Emp_pred = np.array(res)
Emp_pred = cnn.predict(Emp_pred)
print(Emp_pred)

np.savetxt("Emp_ParametersPredictions.txt", Emp_pred)