import numpy as np
import pandas as pd
import os

def mergeRegions(dir = "mergeFiles"):
    os.chdir(dir)
    list_files = []

    pwd = os.getcwd()
    for root, dirs, files in os.walk(pwd):
        for file in files:
            if "csv" in file:
                list_files.append(file)
    print(list_files)
    df = pd.DataFrame()
    for name in list_files:
        df = df.append(pd.read_csv(name, sep = ','))
    df.to_csv("MergedResults.csv", header = True, index = False)





def mergeCal(dir = "mergeFiles"):

    os.chdir("mergeFiles")
    list_files = []

    pwd = os.getcwd()
    for root, dirs, files, in os.walk(pwd):
        for file in files:
            if ".csv" in file and "List" not in file and "results.csv" not in file:
                list_files.append(file)
    f = open('ListInputs.csv', "w")
    print(list_files)
    df = pd.DataFrame()
    for file in list_files:
        name = file[:-4]
        f.write(name)
        df_temp = pd.read_csv(file, sep = ',')
        df["E_"+name] = df_temp['CalibratedE']
    df['cluster_ENG_CALIB_TOT'] = df_temp['cluster_ENG_CALIB_TOT']
    df.to_csv("results.csv", header=True, index = False)




mergeCal()
