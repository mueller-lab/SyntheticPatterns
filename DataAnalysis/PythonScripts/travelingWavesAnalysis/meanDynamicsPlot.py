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

def plot_dynamics(posL, chnls, exDates, timeNum):
    # # make colormap
    viridis = cm.get_cmap('viridis', 256)
    cmapV = viridis(np.linspace(0, 1, 5))

    magma = cm.get_cmap('magma', 256)
    cmapM = magma(np.linspace(0, 1, 5))

    # make a dict to assign cmap to ch
    clrDict = {1:cmapV[3], 2:cmapM[3]}
    lblDict = {1:'sfGFP', 2:'mCherry'}

    # make timeAx
    timeAx = [i*tPerd for i in range(timeNum)] # time in hours
    timeAx_min = 60*np.array(timeAx)

    # make figure
    fig, ax  = plt.subplots(nrows = 1, ncols = 1, figsize= (2.4,1.6))

    chDict = {}
    for ch in chnls:
        yPosAll = []
        for pos in posL:
            for exptDate in exDates:
            # read data and make pandas df
                fName = f"XY{pos}_CH{ch}_measure.csv"
                fPath = inPath/exptDate/'entLawn'/fName
                if fPath.is_file():
                    data_df = pd.read_csv(str(fPath),header=0, index_col = 0)

                    y =data_df['Mean']
                    y = y[:timeNum] # set the last timeNum
                    y = (y-np.min(y))/(np.max(y)-np.min(y))
                    yPosAll.append(y)
        chDict[ch] = np.array(yPosAll)

    for ch in chnls:
        y= chDict[ch]
        yAvg = np.mean(y, axis = 0)
        err = np.std(y, axis=0)
        ax.plot(timeAx, yAvg, marker= ".", markersize=1, color =clrDict[ch], label = lblDict[ch], linewidth = 1)

        # markers, caps, bars = ax.errorbar(timeAx,yAvg, yerr= err, fmt="none", capsize = 0.5, \
        # capthick= 1.0, c = clrDict[ch], elinewidth = 1.)
        #
        # [bar.set_alpha(0.2) for bar in bars] # to make the errorbars transparent
        # [cap.set_alpha(0.2) for cap in caps]

        plt.fill_between(timeAx, yAvg-err, yAvg+err, alpha = 0.2, color = clrDict[ch], edgecolor= 'face')

    ax.grid(False)
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    ax.set_xlabel("Time (h)")
    ax.set_ylabel("Mean intensity (A.U.)")
    # ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    ax.set_ylim([0, 1.06])
    ax.set_xlim(timeAx[0], timeAx[-1])

    ax.tick_params(axis='both', which='major', labelsize=6)
    # ax.set_aspect('equal')

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
    timeNum = 97

    chnls = [1,2]
    posL=['12', '16']
    exDates = ['20230719', '20230720', '20230721']

    plot_dynamics(posL, chnls, exDates, timeNum)

    print("analysis is complete.")
