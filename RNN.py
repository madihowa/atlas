from keras.callbacks import ModelCheckpoint, TensorBoard
from keras.models import Sequential
from keras.layers import Dense, Activation, Flatten, LSTM, SimpleRNN, Embedding
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error
import seaborn as sb
import pandas as pd
import numpy as np
from Graphics import plt_result
import warnings
warnings.filterwarnings('ignore')
warnings.filterwarnings('ignore', category=DeprecationWarning)
warnings.simplefilter(action='ignore', category=FutureWarning)
from xgboost import XGBRegressor
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from CSV_read_write import *
import os
import matplotlib.pyplot as plt
from Normalize import *
import tensorflow as tf
from tensorflow import keras
from keras.utils.vis_utils import plot_model
try:
    # pydot-ng is a fork of pydot that is better maintained.
    import pydot_ng as pydot
except ImportError:
    # pydotplus is an improved version of pydot
    try:
        import pydotplus as pydot
    except ImportError:
        # Fall back on pydot if necessary.
        try:
            import pydot
        except ImportError:
            pydot = None
import math

debug_NN_predict = True
debug_NN_input = False




def FitRNetwork(dir):

    tensorboard_callbacks = TensorBoard(log_dir=dir+"/logs", histogram_freq=1)

    df_test, df_train = Normalize_all()
    df_train = df_train.sample(frac=1) #Mix sample so not all EM or Had are together.

    target = (df_train['cluster_ENG_CALIB_TOT'].to_numpy())

    df_train.drop(['cluster_ENG_CALIB_TOT'],axis = 1 , inplace = True)


#################################################################################################################################
##Cut data down to right observables
##
#################################################################################################################################
    list_input = ReadInputVaribles() #returns currently used observables
    df_train = df_train[list_input]
    df_test = df_test[list_input]

    if debug_NN_input:
        for col in df_train.columns:
            print(col, end = ': ')
            print(df_test.iloc[0][col])



    NN_model = Make_network(df_train.shape[1])


    checkpoint_name = dir+'/Weights-{epoch:03d}--{val_loss:.5f}.hdf5'
    checkpoint = ModelCheckpoint(checkpoint_name, monitor='val_loss', verbose = 1, save_best_only = True, mode ='auto')
    # lr_schedule = keras.callbacks.LearningRateScheduler(learning_schedule)
    callbacks_list = [checkpoint, tensorboard_callbacks]

    #Fit network
    history_callback = NN_model.fit(x=df_train, y=target, shuffle=True, epochs=100, batch_size=128, validation_split = 0.1, callbacks=callbacks_list)
    file = open("NetworkHistory.txt", "a")

    CSV_Callbacks(history_callback)

def NetworkRPredict(dir):
    df1, df2 = Normalize_all()
    if debug_NN_predict == True:
        for col in df1.columns:
            print(col)
    df_test = df1.copy(deep='all')
    df_train = df2


    target = df_train['cluster_ENG_CALIB_TOT']
    df_train.drop(['cluster_ENG_CALIB_TOT'],axis = 1 , inplace = True)


    #############################################################################################################################
    ##Build the Network
    #############################################################################################################################
    NN_model = Make_network(df_train.shape[1])


    #############################################################################################################################
    ##Load the best weights from training
    #############################################################################################################################
    file = find_min_weights(dir)
    min_loss = None
    pwd = os.getcwd()


    wights_file = file # choose the best checkpoint
    NN_model.load_weights(wights_file) # load it

    NN_model = compile_NN(NN_model)

    ##Make predictions
    # predictions = np.exp(NN_model.predict(df_test))
    predictions = (NN_model.predict(df_train))



    #Plot the results
    plt_result(df_train, predictions, target)
    #Save Results as CSV file
    # plt.show()
    plt.savefig(dir + "/figures")



    df_train['cluster_ENG_CALIB_TOT'] = target
    df1, df2 = split()

    df2['CalibratedE'] = predictions
    write_csv_file(df2, dir)

    if debug_NN_predict == True:
        for col in df1.columns:
            print(col)
    tmp = NN_model.layers[0].get_weights()


def FitRNetworkStage2(dir):

    df = pd.read_csv("mergeFiles/results.csv")



    target = df['cluster_ENG_CALIB_TOT']

    df.drop(['cluster_ENG_CALIB_TOT'],axis = 1 , inplace = True)

    NN_model = Make_network(df.shape[1])
    checkpoint_name = dir+'/Weights-{epoch:03d}--{val_loss:.5f}.hdf5'
    checkpoint = ModelCheckpoint(checkpoint_name, monitor='val_loss', verbose = 1, save_best_only = True, mode ='auto')
    callbacks_list = [checkpoint]

    history_callback = NN_model.fit(x=df, y=target, shuffle=True, epochs=10, batch_size=128, validation_split = 0.1, callbacks=callbacks_list)
    file = open("NetworkHistory.txt", "a")

    CSV_Callbacks(history_callback)


def NetworkRPredictStage2(dir):
    df = pd.read_csv("mergeFiles/results.csv")


    target = df['cluster_ENG_CALIB_TOT']
    df.drop(['cluster_ENG_CALIB_TOT'],axis = 1 , inplace = True)


    #############################################################################################################################
    ##Build the Network
    #############################################################################################################################
    NN_model = Make_network(df.shape[1])


    #############################################################################################################################
    ##Load the best weights from training
    #############################################################################################################################
    file = find_min_weights(dir)
    min_loss = None
    pwd = os.getcwd()


    wights_file = file # choose the best checkpoint
    NN_model.load_weights(wights_file) # load it

    NN_model = compile_NN(NN_model)

    ##Make predictions
    # predictions = np.exp(NN_model.predict(df_test))
    predictions = (NN_model.predict(df))



    #Plot the results
    # plt_result(df, predictions, target)
    #Save Results as CSV file
    # plt.show()
    # plt.savefig(dir + "/figures")



    df['cluster_ENG_CALIB_TOT'] = target
    df1, df2 = split()

    df1['CalibratedE'] = predictions
    write_csv_file(df1, dir)


    tmp = NN_model.layers[0].get_weights()


##Here is where you pick the shape and size of the network. This is really where the magic happens. Everything else was just to make runing this work faster.
def Make_network(in_shape):
    NN_model = Sequential()
    # NN_model.add(Embedding(input_dim = in_shape, output_dim = 1, trainable=True))
    # NN_model.add(SimpleRNN(256, activation="relu"))
    NN_model.add(Dense(1024, kernel_initializer='RandomNormal',input_dim = in_shape, activation="relu"))
    # NN_model.add(Dense(2, kernel_initializer='RandomNormal', activation="relu"))
    # NN_model.add(Dense(2, kernel_initializer='RandomNormal', activation="relu"))
    # # NN_model.add(Dense(512, kernel_initializer='RandomNormal', activation="relu"))
    # # NN_model.add(Dense(256, kernel_initializer='RandomNormal', activation="linear"))
    # NN_model.add(Dense(256, kernel_initializer='RandomNormal', activation="linear"))
    NN_model.add(Dense(1024, kernel_initializer='RandomNormal', activation="relu"))
    # NN_model.add(Dense(1024, kernel_initializer='RandomNormal', activation="relu"))
    # NN_model.add(Dense(256, kernel_initializer='RandomNormal', activation="relu"))
    # NN_model.add(Dense(256, kernel_initializer='RandomNormal', activation="linear"))
    # NN_model.add(Dense(256, kernel_initializer='RandomNormal', activation="linear"))
    # NN_model.add(Dense(256, kernel_initializer='RandomNormal', activation="linear"))
    # NN_model.add(Dense(256, kernel_initializer='RandomNormal', activation="linear"))
    # NN_model.add(Dense(256, kernel_initializer='RandomNormal', activation="linear"))
    # NN_model.add(Dense(256, kernel_initializer='RandomNormal', activation="linear"))
    # NN_model.add(Dense(256, kernel_initializer='RandomNormal', activation="linear"))
    # NN_model.add(Dense(256, kernel_initializer='RandomNormal', activation="linear"))
    # NN_model.add(Dense(256, kernel_initializer='RandomNormal', activation="linear"))
    # NN_model.add(Dense(256, kernel_initializer='RandomNormal', activation="linear"))
    # NN_model.add(Dense(256, kernel_initializer='RandomNormal', activation='relu'))
    # NN_model.add(Dense(256, kernel_initializer='RandomNormal', activation='relu'))
    # NN_model.add(Dense(256, kernel_initializer='RandomNormal', activation='relu'))
    # NN_model.add(Dense(256, kernel_initializer='RandomNormal', activation='relu'))
    # NN_model.add(Dense(256, kernel_initializer='RandomNormal', activation='relu'))
    # NN_model.add(Dense(25, kernel_initializer='RandomNormal', activation='relu'))
    # NN_model.add(Dense(25, kernel_initializer='RandomNormal', activation='relu'))
    #



    #
    # for i in range(100):
    #     NN_model.add(Dense(512, kernel_initializer='RandomNormal', activation='relu'))
    #     i = i + 1
    # NN_model.add(Dense(1024, kernel_initializer='normal', activation='sigmoid'))
    NN_model.add(Dense(1, kernel_initializer='RandomNormal',activation=None))
    NN_model = compile_NN(NN_model)
    # NN_model.compile(loss='mean_absolute_error', optimizer='adam', metrics=['mean_absolute_error'])
    NN_model.summary()


    return NN_model

def learning_schedule(epoch):
    learning_rate = .001 * math.exp(-epoch/20)

    tf.summary.scalar('learning rate', data=learning_rate, step=epoch)
    return learning_rate

def compile_NN(NN_model):
    lr_schedule = tf.keras.optimizers.schedules.ExponentialDecay(.001, decay_rate=.36, decay_steps=1e5)
    optimizer = tf.keras.optimizers.Adam(learning_rate=lr_schedule)

    NN_model.compile(loss=tf.keras.losses.MeanAbsolutePercentageError(), optimizer=optimizer, metrics=['mse', 'mae', 'mape'])
    # NN_model.compile(loss=tf.keras.losses.MeanSquaredError(), optimizer=optimizer, metrics=['mse', 'mae', 'mape'])
    return NN_model

##This is not something that helped ever. So that is why it is now just a blank function if you want to do it.
def customActivation(Tensor):
    # tempTensor = Tensor*Tensor*Tensor
    # tempTensor = - tf.exp((-100) * (Tensor * Tensor)) + Tensor
    return tempTensor

def split():
    data_test, data_train  = read_csv_file()

    # data_test = data_test[data_test['cluster_ENG_CALIB_TOT'] >= .5]
    # data_train = data_train[data_train['cluster_ENG_CALIB_TOT'] >= .5]
    Varibles = ReadInputVaribles()
    data_test = data_test[data_test['cluster_SIGNIFICANCE'] >= 2. ]
    data_train = data_train[data_train['cluster_SIGNIFICANCE'] >= 2. ]

    return data_test, data_train
