import pandas as pd
import sys
import csv
import numpy as np
from sklearn.utils import shuffle
import os
from datetime import date
import time
from shutil import copyfile

## List of all the funstions.
__all__ = ["read_csv_file", "find_path", "get_cols_with_no_nans", "find_min_weights", "get_directory", "CSV_Callbacks", "write_csv_file", "ReadInputVaribles", "read_test"]

debug_read = False





def read_csv_file():

    PATH_train, PATH_test = find_path(".csv")
    pd.set_option('display.max_columns', None) # displays all columns
    pd.set_option('display.max_rows', None) # displays all rows

    # drop_columns = ['truthPt', 'truthEta', 'truthPDG', 'cluster_ENG_CALIB_TOT', 'cluster_ENG_CALIB_OUT_T', 'cluster_ENG_CALIB_DEAD_TOT']
    full_train_precut = pd.read_csv(PATH_train, sep = ',')
    # full_train = full_train.drop(columns=drop_columns)
    print("Done reading train")
    full_test_precut = pd.read_csv(PATH_test, sep = ',')
    # full_test = full_test.drop(columns=drop_columns)
    print("Done reading test")
    full_train = full_train_precut
    full_test = full_test_precut
    # full_train = full_train_precut[full_train_precut['cluster_SIGNIFICANCE'] >= 0.]
    # full_test = full_test_precut[full_test_precut['cluster_SIGNIFICANCE'] >= 0.]

        # df2 = df2[df2['cluster_ENG_CALIB_TOT'] >= 10]
        # df1 = df1[df1['cluster_ENG_CALIB_TOT'] >= 10]
    ##add column for use in testing

    ##Print column names and length of data if debug is set to true.
    if debug_read:
        print("Total entries in file: ")
        #print("Total entries in file: ", end=='')
        print(len(full_train))
        for col in full_test.columns:
            print(col)
            ##print(col, end=='\t')
        print("")
        for col in full_test.columns:
            ##print((full_test.iloc[0][col]), end=='\t')
            print((full_test.iloc[0][col]))
        print("")
    return full_test, full_train

def read_test():
    PATH_train, PATH_test = find_path(".csv")
    drop_columns = ['truthPt', 'truthEta', 'truthPDG', 'cluster_ENG_CALIB_TOT', 'cluster_ENG_CALIB_OUT_T', 'cluster_ENG_CALIB_DEAD_TOT']
    full_test = pd.read_csv(PATH_test, sep = ',')
    full_test = full_test.drop(columns=drop_columns)
    print("Done reading test")

    return full_test

##Return the path of the test and train csv files. Note that it finds all CSV files with train and .csv in them and having more then one will cause issues same with test and .csv
def find_path(str_find):
    pwd = os.getcwd()
    for root, dirs, files in os.walk(pwd):
        for file in files:
            if "EM" in root:
                if file.endswith(str_find):
                    if "train" in file or "Train" in file:
                        if debug_read:
                            print(file)
                        file_dir_train = os.path.join(root, file)
                    if "test" in file or "Test" in file:
                        if debug_read:
                            print(file)
                        file_dir_test = os.path.join(root, file)

    return file_dir_train, file_dir_test


def get_cols_with_no_nans(df,col_type):
    '''
    Arguments :
    df : The dataframe to process
    col_type :
          num : to only get numerical columns with no nans
          no_num : to only get nun-numerical columns with no nans
          all : to get any columns with no nans
    '''
    if (col_type == 'num'):
        predictors = df.select_dtypes(exclude=['object'])
    elif (col_type == 'no_num'):
        predictors = df.select_dtypes(include=['object'])
    elif (col_type == 'all'):
        predictors = df
    else :
        print('Error : choose a type (num, no_num, all)')
        return 0
    cols_with_no_nans = []
    for col in predictors.columns:
        if not df[col].isnull().any():
            cols_with_no_nans.append(col)
    return cols_with_no_nans


##This file returns the set of weights in a specific run that proformed the best. That is if over fitting occurs then it will pick the set before over fitting occurs. However this really never happened to me.
def find_min_weights(directory):
        pwd = os.getcwd()
        min_loss = None
        for root, dirs, files in os.walk(directory):
            for file in files:
                if file.endswith('hdf5'):
                    splt = file.split('--')
                    if min_loss == None:
                            min_loss=float(splt[1][:-5])
                            min_file = os.path.join(root, file)
                    elif min_loss > float(splt[1][:-5]):
                            min_loss = float(splt[1][:-5])
                            min_file = os.path.join(root, file)
                    # print(file)
                    # print(splt[1][:-5])
                    # print(float(splt[1][:-5]))
        return min_file


def get_directory(day, run):
    """
    Inputs: day->string , run->int
    Output: directory->string
    Summary: creates a directory with stylized name and copies the ListInputs.csv into this directory
    """
    if day == "today":
        day = (date.today()).isoformat()
    else:
        try:
            date.fromisoformat(day)
        except:
            print("error please input date as YYYY-MM-DD")
            quit()

    pwd = os.getcwd()
    print(pwd)
    runs = 1
    dir_exists = False
    print("Run number is: ", run)
    if run == -1:
        for root, dirs, files in os.walk(pwd):
            for dir in dirs:
                if day in dir:
                    runs = runs + 1
        directory = "Results_"+day+"_run_"+str(runs)
        os.mkdir(os.path.join(pwd, directory))
        print(directory)
        copyfile("ListInputs.csv", os.path.join(directory,"ListInputs.csv"))
    else:
        directory = "Results_"+day+"_run_"+str(run)
        for root, dirs, files in os.walk(pwd):
            for dir in dirs:
                if directory in dir:
                    dir_exists=True
                    break
        if not dir_exists:
            print("error that run doesn't exist")
            quit()
    return directory

def write_csv_file(df, dir, file_name='results.csv'):
    print("testing the write func:"  + dir)
    df.to_csv(os.path.join(dir, file_name), header=True, index=False)
    print("CSV file written")
    return

def ReadInputVaribles():
    list_input = []

    file = 'ListInputs.csv'


    with open(file) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter = ',')
            for row in csv_reader:
                if '#' in row[0]: continue
                list_input = list_input + row
    # print(list_input)
    return list_input

def CSV_Callbacks(callbacks):

    f = open("callback_history.csv", "w")
    print(type(callbacks.history))
    for title in callbacks.history:
        print(title)
        ##print(title, end==',')

    history_df = pd.DataFrame.from_dict(callbacks.history)
    history_df.to_csv("callback_history.csv", index=False)
