import pandas as pd
import numpy as np
from CSV_read_write import *
__Debug__ = True
__all__ = ["Normalize_all", "small_set"]


##This function was used to try many different ways of normalizing the data to try and get better results.
def Normalize_all():
    from RNN import split


    df1, df2 = split()
    list_input = ReadInputVaribles()
    # list_input.append("EM_Shower")

    # df2 = df2_pre[df2_pre['cluster_ENG_CALIB_TOT'] >= 0]
    # df1 = df1_pre[df1_pre['cluster_ENG_CALIB_TOT'] >= 0]
    # df1['EM_Pro'] = np.exp(df1['EM_Pro'])
    # df2['EM_Pro'] = np.exp(df2['EM_Pro'])
    #
    # df1['EM_Pro'] = np.exp(df1['EM_Pro'])
    # df2['EM_Pro'] = np.exp(df2['EM_Pro'])
    #
    # df1['cluster_CENTER_LAMBDA'] = np.square(df1['cluster_CENTER_LAMBDA'])
    # df2['cluster_CENTER_LAMBDA'] = np.square(df2['cluster_CENTER_LAMBDA'])

    if __Debug__ == True:
        for col in df1.columns:
            print(col)

    # df1['clusterE'] = np.log(df1['clusterE'])
    # df2['clusterE'] = np.log(df2['clusterE'])
    df_test = df1[list_input]
    df_train = df2[list_input]

    # n_col = len(df_test.columns)
    # list_mean_EM = []
    # list_mean_Had = []
    # list_col_names = []
    # print("        ")
    for col in df_test.columns:
        print(col)


    # print(df_em.mean())
    # print(df_em.std())

    # print(df_had.mean())
    # print(df_had.std())
    # print("Cluster Min is: ", df_train['cluster_CENTER_LAMBDA'].min())
    # print(df_test[['cluster_FIRST_ENG_DENS', 'cluster_CENTER_LAMBDA']][:10])
    # print("Clusters with zero lambda: ", (df_train.loc[df_train['cluster_CENTER_LAMBDA'] == 0]).count())
    # df_test = df_test[df_test['cluster_CENTER_LAMBDA'] != 0]
    # df_test = df_test[df_test['cluster_FIRST_ENG_DENS'] != 0]
    # df_train = df_train[df_train['cluster_CENTER_LAMBDA'] != 0]
    # df_train = df_train[df_train['cluster_FIRST_ENG_DENS'] != 0]
    # if "cluster_FIRST_ENG_DENS" in list_input:
    #      df_test['cluster_FIRST_ENG_DENS'] = np.log10(df_test['cluster_FIRST_ENG_DENS'])
    #      print("Taking Log of first energy density")
    #      df_train['cluster_FIRST_ENG_DENS'] = np.log10(df_train['cluster_FIRST_ENG_DENS'])
    # if 'cluster_CENTER_LAMBDA' in list_input:
    #      df_test['cluster_CENTER_LAMBDA'] = np.log10(df_test['cluster_CENTER_LAMBDA'])
    #      df_train['cluster_CENTER_LAMBDA'] = np.log10(df_train['cluster_CENTER_LAMBDA'])
    #      print("Taking Log of lambda")
    # print(df_test[['cluster_FIRST_ENG_DENS', 'cluster_CENTER_LAMBDA']][:10])




    df_test_norm = df_test.subtract(df_test.mean())
    df_test_norm = df_test_norm.divide(df_test.std())
    # df_test_norm = df_test.subtract(0)
    # df_test_norm = df_test_norm.divide(1)

    df_train_norm = df_train.subtract(df_train.mean())
    df_train_norm = df_train_norm.divide(df_train.std())
    # df_train_norm = df_train.subtract(0)
    # df_train_norm = df_train_norm.divide(1)

    # df_test_norm['cluster_ENG_CALIB_TOT'] = np.log(df1['cluster_ENG_CALIB_TOT'])
    # df_train_norm['cluster_ENG_CALIB_TOT'] = np.log(df2['cluster_ENG_CALIB_TOT'])

    # df_test_norm['clusterE'] = (df1['clusterE'])
    # df_train_norm['clusterE'] = (df2['clusterE'])
    df_train_norm['cluster_ENG_CALIB_TOT'] = df2['cluster_ENG_CALIB_TOT']
    df_test_norm['cluster_ENG_CALIB_TOT'] = df1['cluster_ENG_CALIB_TOT']

    #
    # df_test_norm['EM_Shower'] = df1['EM_Shower']
    # df_test_norm['cluster_EM_PROBABILITY'] = df1['cluster_EM_PROBABILITY']
    # df_train_norm['EM_Shower'] = df2['EM_Shower']
    # df_train_norm['cluster_EM_PROBABILITY'] = df2['cluster_EM_PROBABILITY']

    #######################################################
    #Save Data for seperated Had and EM.
    # df_em = df_train[df_train['EM_Shower'] == 1]
    # df_had = df_train[df_train['EM_Shower'] == 0]
    # write_csv_file(df_em, '', "EM_showers.csv")
    # write_csv_file(df_had, '', "Had_showers.csv")

    print("The mean of the test is \n", df_test.mean())
    print("The std of the test is \n", df_test.std())
    print("The mean of the training is \n", df_train.mean())
    print("The std of the training is \n", df_train.std())

    return df_test_norm, df_train_norm

# Normalize_all()


##This makes a data set with only 10 values so that other code could be tested quickly.
def small_set():
    from RNN import split

    df_test, df_train = Normalize_all()
    n = 10
    df_test = df_test[:n]
    df_train = df_train[:n]
    print(df_train)
    return df_test, df_train
