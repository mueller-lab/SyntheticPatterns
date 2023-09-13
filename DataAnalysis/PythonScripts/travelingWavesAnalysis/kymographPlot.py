"""
Amit Landge
Custom code to analyse and plot data from the plate reader assays.

"""
import os
import time
from datetime import datetime
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from scipy.optimize import curve_fit
# import argparse

# functions
def smoothConvolve(data, sz = 5, mode = "same"):
    kernel = np.ones(sz)/sz
    smthData = np.convolve(data, kernel, mode= mode)
    return smthData

def plot_kymo(xypos, ch, timeNum):
    # read data and make pandas df
    fName = f"radProf_XY{xypos}_CH{ch}.csv"
    fPath = inPath/fName
    data_df = pd.read_csv(str(fPath),header=0, index_col = 0)

    # exclude the edge of the lawn (eg. radius > 5)
    data_df = data_df.loc[data_df.index<5.0, :]

    # process the data and make space-time plots
    timeAx = [i*tPerd for i in range(timeNum)] # time in hours

    data_df.columns = timeAx
    lenAx = data_df.index

    # print(data_df.head) # dataframe looks good

    # vmax = np.amax(data_df.values)
    # vmin = np.amin(data_df.values)

    # make numpy array of data_df
    dataArr0 = np.array(data_df)
    rows, cols = dataArr0.shape

    # # make the data smooth
    # valsSmth = np.zeros((rows, cols))
    # for j in range(cols):
    #     jCol = smoothConvolve(vals[:,j])
    #     valsSmth[:,j]=jCol

    #repeat each time axis 10-fold to make the space-time plot visualization better
    tMulti = 2
    dataArr = np.zeros((rows, tMulti*cols))
    for j in range(cols):
        for k in range(tMulti):
            dataArr[:,tMulti*j+k] = dataArr0[:,j]

    # make a dict to assign cmap to ch
    cmapDict = {1:'viridis', 2:'magma', 4:'gray'}

    # make figure
    fig, ax  = plt.subplots(nrows = 1, ncols = 1, figsize= (2.0,2.0))
    pos = ax.imshow(dataArr, origin = 'lower', cmap =cmapDict[ch]) #cmap = 'magma'

    ax.grid(False)
    ax.set_ylabel("r (mm)")
    ax.set_xlabel("Time (h)")

    lenSt = 2
    timeSt = 10

    # set x and y-ticks
    yTic = [int(lenSt*i/lenAx[0]) for i in range(1+int(lenAx[-1]/lenSt))]
    y_labels = [lenSt*i for i in range(1+int(lenAx[-1]/lenSt))]
    ax.set_yticks(yTic)
    ax.set_yticklabels(y_labels)

    xTic = [int(timeSt*i*tMulti/timeAx[1]) for i in range(1+int(timeAx[-1]/timeSt))]
    x_labels = [timeSt*i for i in range(1+int(timeAx[-1]/timeSt))]
    ax.set_xticks(xTic)
    ax.set_xticklabels(x_labels)

    # ax.set_ylim(timeAx[0], timeAx[-1])
    # ax.set_xlim(xVals[0], xVals[-1])
    ax.tick_params(axis='both', which='major', labelsize=6)

    # add colorbar
    plt.colorbar(pos, fraction = 0.04, pad = 0.3, orientation = 'horizontal') # label = 'Intensity val.',

    plt.tight_layout()

    #remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    figName_0 = f"XY{xypos}_CH{ch}_SpaceTime.png"
    figName_1 = f"XY{xypos}_CH{ch}_SpaceTime.pdf"

    # plt.show()
    plt.savefig(outPath/figName_0, transparent=True, dpi = 300)
    plt.savefig(outPath/figName_1, transparent=True)
    plt.close()

# main code
if __name__=="__main__":
    print("starting analysis ...")

    # global plot settings
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["font.size"] = 6
    plt.rcParams['legend.fontsize'] = 6

    cwdPath=Path(os.path.abspath(os.getcwd()))
    inPath= cwdPath
    datetime_str=(datetime.now().strftime("%Y%m%d_"))
    outDir = str(datetime_str+'kymographs')
    outPath=inPath/outDir
    outPath.mkdir(parents=True, exist_ok=True)

    tPerd = 0.25
    timeNum = 105

    for ch in [1,2,4]:
        for xypos in [17,18]:
            plot_kymo(xypos, ch, timeNum)

    print("analysis is complete.")
