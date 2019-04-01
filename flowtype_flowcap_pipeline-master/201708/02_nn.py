# aya43@sfu.ca 20161223

import numpy as np
import pandas as pd
import os
import model
import keras.backend as K
from keras.optimizers import SGD
from keras.optimizers import Adadelta
from utils import LearnRateScheduler, WeightsWriter
from keras.callbacks import EarlyStopping
import cPickle as pickle
import pdb

#Input
results_dir = './results'
matrix = pd.read_csv(results_dir + '/matrixCountAdj.csv', sep=',', header=True)
#matrix = pd.read_csv('./results/matrixCount.csv', sep=',', header=True)
#matrix = pd.read_csv('./results/matrixProp.csv', sep=',', header=True)
sampleMeta = pd.read_csv('./results/sampleMeta_features.csv', sep=',', header=True)

#Output
results_dir = results_dir + '/keras_dual'
if not os.path.exists(results_dir):
    os.makedirs(results_dir)
best_weights_file = results_dir + '/best_weights.pkl'
final_weights_file = results_dir + '/final_weights.pkl'
history_file = results_dir + '/history.pkl'

#Training
num_categories = 2
do = .5 #dropout rate
learning_rate = .001
momentum = 0.9
early_stop_patience = 3
lr_step_decay = 0.8
lr_decay_patience = 2
num_epcs = 25

#Prepare Input
#subset sampleMeta columns

#Compile Model
model = model.dual(pheno_size, num_categories, do)
ad = Adadelta()
sgd = SGD(lr=learning_rate, momentum=momentum, nesterov=True, clipnorm=5.)
model.compile(loss='categorical_crossentropy', optimizer=ad, metrics=['accuracy'])

wwriter = WeightsWriter(filepath=best_weights_file, monitor='val_acc', verbose=0, save_best_only=True, mode='auto')
early_stop = EarlyStopping(monitor='val_acc', patience=early_stop_patience, verbose=0, mode='auto')
lr_schd = LearnRateScheduler(lr_step_decay, monitor='val_acc', patience=lr_decay_patience, verbose=0, mode='auto')

hist = model.fit(
        [], target_train,
        samples_per_epoch=2000,
        nb_epoch=num_epcs,
        validation_data = [[],target_val],
        nb_val_samples=1000,
        callbacks=[early_stop, lr_schd, wwriter])

weights = model.get_weights()
with open(final_weights_file, 'wb') as f:
    pickle.dump(weights, f)

try:
    hist.history['lr']=K.get_value(hist.model.optimizer.lr)
    with open(history_file, 'w') as f:
        pickle.dump(hist.history, f)
except Exception as e:
    print e
    pdb.set_trace()

