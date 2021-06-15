from CSV_read_write import *
import glob, os, argparse, codecs, sys
import pandas as pd
from RNN import NetworkRPredictStage2, FitRNetworkStage2
import matplotlib.pyplot as plt
import shutil, os

#Set the right input files
parser = argparse.ArgumentParser(description='Creates a NN and trains it to calibrate topo-cluster energies')
parser.add_argument("-r", type=int, help="Pick a run number to call previous data else will make new network.", default=-1)
parser.add_argument("-d", type=str, help="Pick a date which the run is from, default today", default="today")
args = parser.parse_args()
run = args.r
day = args.d

if run == -1:
    new_run = True
else:
    new_run = False

directory = get_directory(day, run)
file = open("NetworkHistory.txt", "a")
if run == -1:
    run = 0
    for root, dirs, files in os.walk(os.getcwd()):
        for dir in dirs:
            if day in dir:
                run = run + 1
file.write('\n{:}\t|{:1}\t|'.format(day, run))

shutil.copy("RNN.py", directory)
new_run = True
if new_run:
    FitRNetworkStage2(directory)
else:
    file.write('\t\t')

NetworkRPredictStage2(directory)
