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
import sklearn.metrics as metrics
from sklearn.metrics import log_loss
from sklearn.metrics import confusion_matrix

batch_size = 500

epochs = 500

num_classes = 3

def create_cnn(xtest):
	inputShape = (xtest.shape[1], xtest.shape[2])
	inputs = Input(shape=inputShape)
	x = inputs
	x = Conv1D(125, kernel_size=3, activation='relu', use_bias=False, input_shape=(xtest.shape[1], xtest.shape[2]))(x)
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
	# The final fully-connected layer head will have a softmax dense layer
	logits = Dense(num_classes, name='logits')(x)
	out = Activation('softmax')(logits)

	# Construct the CNN
	model = Model(inputs, out)
	# Return the CNN
	return model


def fit_TemperatureCalibration(train_X_y, valid_X_y=None, epochs=100):
    ### From: https://github.com/cerlymarco/MEDIUM_NoteBook/blob/master/NeuralNet_Calibration/NeuralNet_Calibration.ipynb;
    ###inspired by: https://github.com/stellargraph/stellargraph/blob/develop/stellargraph/calibration.py ###
    
    T = tf.Variable(tf.ones(shape=(1,)), name="T")
    history = []
    early_stopping = False
    optimizer = SGD(learning_rate=0.001)
    
    def cost(T, x, y):

        scaled_logits = tf.multiply(x=x, y=1.0 / T)

        cost_value = tf.reduce_mean(
            tf.nn.softmax_cross_entropy_with_logits(logits=scaled_logits, labels=y)
        )

        return cost_value

    def grad(T, x, y):

        with tf.GradientTape() as tape:
            cost_value = cost(T, x, y)

        return cost_value, tape.gradient(cost_value, T)
    
    
    X_train, y_train = train_X_y
    if valid_X_y:
        X_valid, y_valid = valid_X_y
        early_stopping = True
    
    
    for epoch in range(epochs):
        train_cost, grads = grad(T, X_train, y_train)
        optimizer.apply_gradients(zip([grads], [T]))
        if early_stopping:
            val_cost = cost(T, X_valid, y_valid)
            if (len(history) > 0) and (val_cost > history[-1][1]):
                break
            else: 
                history.append([train_cost, val_cost, T.numpy()[0]])
        else:
            history.append([train_cost, T.numpy()[0]])

    history = np.asarray(history)
    temperature = history[-1, -1]
    
    return temperature

def cal_softmax(x):
    """
    From: https://github.com/markus93/NN_calibration
    Compute softmax values for each sets of scores in x.
    
    Parameters:
        x (numpy.ndarray): array containing m samples with n-dimensions (m,n)
    Returns:
        x_softmax (numpy.ndarray) softmaxed values for initial (m,n) array
    """
    e_x = np.exp(x - np.max(x))
    return e_x / e_x.sum(axis=1, keepdims=1)



def evaluate(probs, y_true, verbose = False, normalize = False, bins = 15):
    """
    Modified from: https://github.com/markus93/NN_calibration
    Evaluate model using various scoring measures: Error Rate, ECE, MCE, NLL, Brier Score
    
    Params:
        probs: a list containing probabilities for all the classes with a shape of (samples, classes)
        y_true: a list containing the actual class labels
        verbose: (bool) are the scores printed out. (default = False)
        normalize: (bool) in case of 1-vs-K calibration, the probabilities need to be normalized.
        bins: (int) - into how many bins are probabilities divided (default = 15)
        
    Returns:
        (error, loss, brier), returns various scoring measures
    """
    
    preds = np.argmax(probs, axis=1)  # Take maximum confidence as prediction
    
    if normalize:
        confs = np.max(probs, axis=1)/np.sum(probs, axis=1)
        # Check if everything below or equal to 1?
    else:
        confs = np.max(probs, axis=1)  # Take only maximum confidence
    
    accuracy = metrics.accuracy_score(y_true, preds) * 100
    error = 100 - accuracy
        
    loss = log_loss(y_true=y_true, y_pred=probs)
        
    if verbose:
        print("Accuracy:", accuracy)
        print("Error:", error)
        print("Loss:", loss)
    
    return (error, loss)
    
#SSCont = np.load("../Mod_SSCont.npy",mmap_mode='r')
#SSBackCol = np.load("../Mod_SSBackCol.npy",mmap_mode='r')
SSHBackCol = np.load("../Mod_SSHBackCol.npy",mmap_mode='r')
SSHclim = np.load("../Mod_SSHclim.npy",mmap_mode='r')
#SSHwinds = np.load("../Mod_SSHwinds.npy",mmap_mode='r')
EastWest = np.load("../Mod_EastWest.npy",mmap_mode='r')

x=np.concatenate((SSHBackCol[:,0:1000,:],SSHclim[:,0:1000,:],EastWest[:,0:1000,:]),axis=0)

#transform major alleles in -1 and minor 1
for arr,array in enumerate(x):
  for idx,row in enumerate(array):
    if np.count_nonzero(row) > len(row)/2:
      x[arr][idx][x[arr][idx] == 1] = -1
      x[arr][idx][x[arr][idx] == 0] = 1
    else:
      x[arr][idx][x[arr][idx] == 0] = -1

y=[0 for i in range(len(SSHBackCol))]
y.extend([1 for i in range(len(SSHclim))])
y.extend([2 for i in range(len(EastWest))])
#y.extend([3 for i in range(len(SSHclim))])
#y.extend([4 for i in range(len(EastWest))])

y = np.array(y)

del (SSHBackCol,SSHclim,EastWest)
print (len(x), len(y))
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

x=np.array(x)
#x=np.swapaxes(x,1,2)

xtrain, xtest = x[int(len(y)*.25):], x[:int(len(y)*.25)]
ytrain, ytest = y[int(len(y)*.25):], y[:int(len(y)*.25)]

ytest = np_utils.to_categorical(ytest, num_classes)
ytrain = np_utils.to_categorical(ytrain, num_classes)

# Create the CNN network
model = create_cnn(xtest)

opt = SGD(learning_rate=0.001)

model.compile(loss=keras.losses.categorical_crossentropy,
	              optimizer=opt,
	              metrics=['accuracy'])

model.summary()

earlyStopping = EarlyStopping(monitor='val_accuracy', patience=150, verbose=0, mode='max', restore_best_weights=True)
start = time.time()
model.fit(xtrain, ytrain, batch_size=batch_size,
          epochs=epochs,
          verbose=1,
          validation_data=(xtest, ytest),callbacks=[earlyStopping])
print (f'Time: {time.time() - start}')

model.save(filepath='Trained_uncalibrated_Euphorbia.acc.mod')

#SSCont = np.load("../testSet/Mod_SSCont.npy",mmap_mode='r')
#SSBackCol = np.load("../testSet/Mod_SSBackCol.npy",mmap_mode='r')
SSHBackCol = np.load("../testSet/Mod_SSHBackCol.npy",mmap_mode='r')
SSHclim = np.load("../testSet/Mod_SSHclim.npy",mmap_mode='r')
#SSHwinds = np.load("../testSet/Mod_SSHwinds.npy",mmap_mode='r')
EastWest = np.load("../testSet/Mod_EastWest.npy",mmap_mode='r')

xcal=np.concatenate((SSHBackCol[0:5000,0:1000,],SSHclim[0:5000,0:1000,],EastWest[0:5000,0:1000,]),axis=0)
xpred=np.concatenate((SSHBackCol[5000:10000,0:1000,],SSHclim[5000:10000,0:1000,],EastWest[5000:10000,0:1000,]),axis=0)

#transform major alleles in -1 and minor 1
for arr,array in enumerate(xcal):
  for idx,row in enumerate(array):
    if np.count_nonzero(row) > len(row)/2:
      xcal[arr][idx][xcal[arr][idx] == 1] = -1
      xcal[arr][idx][xcal[arr][idx] == 0] = 1
    else:
      xcal[arr][idx][xcal[arr][idx] == 0] = -1

#transform major alleles in -1 and minor 1
for arr,array in enumerate(xpred):
  for idx,row in enumerate(array):
    if np.count_nonzero(row) > len(row)/2:
      xpred[arr][idx][xpred[arr][idx] == 1] = -1
      xpred[arr][idx][xpred[arr][idx] == 0] = 1
    else:
      xpred[arr][idx][xpred[arr][idx] == 0] = -1


y=[0 for i in range(len(SSHBackCol[0:5000,:,:]))]
y.extend([1 for i in range(len(SSHclim[0:5000,:,:]))])
y.extend([2 for i in range(len(EastWest[0:5000,:,:]))])
#y.extend([3 for i in range(len(SSHclim[0:5000,:,:]))])
#y.extend([4 for i in range(len(EastWest[0:5000,:,:]))])

y = np.array(y)

del (SSHBackCol,SSHclim,EastWest)
print (len(xcal), len(xpred), len(y))
shf = list(range(len(xcal)))
shuffle(shf)

y_shf = y[shf]
xcal = xcal[shf]

missD_perc = 4.4
#Add missing data as 0s, according to a specifies missing data percentage
missD = int(xcal.shape[1]*xcal.shape[2]*(missD_perc/100))
for i in range(xcal.shape[0]):
  for m in range(missD):
    j = random.randint(0, xcal.shape[1] - 1)
    k = random.randint(0, xcal.shape[2] - 1)
    xcal[i][j][k] = 0
del(missD)

#Add missing data as 0s, according to a specifies missing data percentage
missD = int(xpred.shape[1]*xpred.shape[2]*(missD_perc/100))
for i in range(xpred.shape[0]):
  for m in range(missD):
    j = random.randint(0, xpred.shape[1] - 1)
    k = random.randint(0, xpred.shape[2] - 1)
    xpred[i][j][k] = 0
del(missD)

xcal=np.array(xcal)
xpred=np.array(xpred)
#xcal=np.swapaxes(xcal,1,2)
#xpred=np.swapaxes(xpred,1,2)

xtrain, xtest = xcal[int(len(y)*.5):], xcal[:int(len(y)*.5)]
ytrain, ytest = y_shf[int(len(y_shf)*.5):], y_shf[:int(len(y_shf)*.5)]

ytest = np_utils.to_categorical(ytest, num_classes)
ytrain = np_utils.to_categorical(ytrain, num_classes)

model_score = Model(inputs=model.input, outputs=model.layers[-1].output)
X_train_calib = model_score.predict(xtrain)
X_valid_calib = model_score.predict(xtest)
temperature = fit_TemperatureCalibration((X_train_calib,ytrain), (X_valid_calib,ytest), epochs=100)
print ("Temperature",temperature)

pred = model.predict(xpred)
pred_cat = [i.argmax() for i in pred]
print (confusion_matrix(y, pred_cat))
print (confusion_matrix(y, pred_cat) / float(len(y)))

np.savetxt("Pred_testSet_UncalibratedModelPredictions.txt", pred)

pred_cal = model_score.predict(xpred)
scaled_prediction = cal_softmax(pred_cal/temperature)

scpred_cat = [i.argmax() for i in scaled_prediction]
print (confusion_matrix(y, scpred_cat))
print (confusion_matrix(y, scpred_cat) / float(len(y)))

np.savetxt("Pred_testSet_CalibratedModelPredictions.txt", scaled_prediction)


infile=np.loadtxt("../input_Euphorbia.txt")
inp=np.array(infile)
num_samples=100
res = []
for i in range(0,num_samples):
	idx = np.random.choice(inp.shape[0], 1000, replace=False)
	n = inp[idx,:]
	res.append(np.array(n))

Emp = np.array(res)
Emp_pred = model.predict(Emp)

np.savetxt("Pred_Emp_UnCalModelPredictions.txt", Emp_pred)


Emp_calpred = model_score.predict(Emp)
Emp_calpred = cal_softmax(Emp_calpred/temperature)

np.savetxt("Pred_Emp_CalModelPredictions.txt", Emp_calpred)

print("Uncalibrated Error %f; loss %f" % evaluate(pred, y, verbose=False, normalize=True))
print("Calibrated Error %f; loss %f" % evaluate(scaled_prediction, y, verbose=False, normalize=True))