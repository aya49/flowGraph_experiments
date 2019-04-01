# aya43@sfu.ca 20161227

# Create Keras Models

from keras.models import Model, Sequential
from keras.layers import Flatten, Dense, Input, merge, Dropout
from keras.layers import Convolution2D, MaxPooling2D
from keras import backend as K


#vector two-branch model
def dual(pheno_size, meta_size, num_categories, do):
    pheno_inp=Input(shape=pheno_size, name='pheno')
    x=Dense(pheno_size[0], activation='relu')(pheno_inp)
    x=Dropout(do)(x)
    x=Dense(pheno_size[0], activation='relu')(x)
    x=Dropout(do)(x)
    x=Dense(1, activation='sigmoid')(x)
    model=Model(input=inp, output=x)
    return model

#vector single-branch model
def single(pheno_size, num_categories):
    pheno_inp=Input(shape=(3,)+img_size, name='img')
    inp = Input(shape=vec_size)
    x = Dense(4096, activation='relu')(inp)
    x=Dropout(do)(x)
    x=Dense(4096, activation='relu')(x)
    x=Dropout(do)(x)
    x=Dense(1, activation='sigmoid')(x)
    model=Model(input=inp, output=x)
    return model