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

def plot_dynamics(chnls, timeNum):
    # make a dict to assign cmap to ch
    clrDict = {1:'g', 2:'r'}
    lblDict = {1:'sfGFP', 2:'mCherry'}

    # make timeAx
    timeAx = [i*tPerd for i in range(timeNum)] # time in hours

    # make figure
    fig, ax  = plt.subplots(nrows = 1, ncols = 1, figsize= (2.5,1.5))

    for ch in chnls:
        # read data and make pandas df
        fName = f"CH{ch}_measure.csv"
        fPath = inPath/fName
        data_df = pd.read_csv(str(fPath),header=0, index_col = 0)

        y =data_df['Mean']
        y = (y-10.0)/(np.max(y)-10)

        ax.plot(timeAx, y, marker= ".", markersize=1, color =clrDict[ch], label = lblDict[ch], linewidth = 1)

    ax.grid(False)
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    ax.set_xlabel("Time (h)")
    ax.set_ylabel("Mean intensity (A.U.)")
    # ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    ax.set_ylim([0, 1.02])
    ax.set_xlim(timeAx[0], timeAx[-1])

    ax.tick_params(axis='both', which='major', labelsize=6)

    plt.tight_layout()

    #remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    figName_0 = f"mean_vs_time.png"
    figName_1 = f"mean_vs_time.pdf"

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
    outDir = str(datetime_str+'meanDynamics')
    outPath=inPath/outDir
    outPath.mkdir(parents=True, exist_ok=True)

    tPerd = 0.25
    timeNum = 105

    chnls = [1,2]
    plot_dynamics(chnls, timeNum)

    print("analysis is complete.")
