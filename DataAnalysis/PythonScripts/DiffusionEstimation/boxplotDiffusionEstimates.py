
"""
Amit Landge
Custom code to analyse and plot data.

"""
import os
import time
from datetime import datetime
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
# from scipy.stats import linregress
import seaborn as sns
# import argparse

# functions

# main code
if __name__=="__main__":
    print("starting analysis ...")

    # global plot settings
    plt.rcParams["font.family"] = "Arial"
    # plt.rcParams["font.style"] = "normal"
    # plt.rcParams["font.weight"] = 100
    plt.rcParams["font.size"] = 6
    plt.rcParams['legend.fontsize'] = 6
    plt.rcParams['pdf.fonttype'] = 42

    cwdPath=Path(os.path.abspath(os.getcwd()))
    inPath= cwdPath
    datetime_str=(datetime.now().strftime("%Y%m%d_")) # _%H%M%S
    outDir = str(datetime_str+"output_diffEsti_BoxPlot")
    outPath=inPath/outDir
    outPath.mkdir(parents=True, exist_ok=True)

    # structure the data and save as a .csv file
    dataDict={}
    mols = [r'Dextran_{3kDa}', r'Dextran_{10kDa}', r'AHL_6', r'AHL_{12}']

    Dex3_all = [122.17,129.2,121.28,114.29,113.63,116.95,120.93] # n = 7
    Dex10_all = [89.45,66.17,73.71,51.51,58.39,61.15,67.52] # n = 7
    D6_all = [472.00, 425.00, 339.00, 310.00, 276.00, 290.00, 302.00, 275.00, 279.08, 277.03, 245.66, 260.14] # n = 12
    D12_all = [443.37, 392.6, 300.74, 309.33, 283.24, 278.83, 288.94, 276.14, 247.58, 365.45, 329.91] # n = 11

    Ds = [Dex3_all, Dex10_all, D6_all, D12_all]

    # add 'NaN' to make all the lists of same length (12)
    for i, D_all in enumerate(Ds):
        if len(D_all)<12:
            D_all = D_all + (12-len(D_all))*[np.nan]
            Ds[i] = D_all


    for i,m in enumerate(mols):
        dataDict[m] = Ds[i]
        i+=1

    diffu_df = pd.DataFrame(dataDict)
    diffu_df.columns.set_names(['mol'], inplace=True)
    diffu_df = diffu_df.stack(level=0)
    diffu_df.index.set_names(['id', 'mol'], inplace=True)
    diffu_df= diffu_df.to_frame()
    diffu_df.columns = ['diffCoef']

    diffu_df= diffu_df.reset_index()
    diffu_df = diffu_df[['mol', 'diffCoef']]

    diffu_df = diffu_df.set_index(['mol'])
    # diffu_df = diffu_df.sort_index()
    diffu_df= diffu_df.loc[mols,:]
    diffu_df.to_csv(outPath/"diffusionEstimates.csv")

    # get mean and standard deviation
    meanStd_df = diffu_df.groupby('mol', sort= False).agg({'diffCoef':['mean','std']})
    meanStd_df.to_csv(outPath/"meanStd_df.csv")


    # process the data
    # steps - 1) make boxplots
    diffu_df = diffu_df.reset_index()
    diffu_df['diffCoef'] = diffu_df['diffCoef'].astype(float) # conver column to float

    # x-positions for the plot
    xScat = 7*[1]+7*[3]+12*[5]+11*[7]
    # xScat = [i for j in xScat0 for i in j]
    diffu_df['xPos']=xScat

    #plot fold-change as box-plot with scatter
    fig, ax  = plt.subplots(nrows = 1, ncols = 1, figsize= (1.6,1.6))

    xlabs = [r'$Dex_{3kDa}$', r'$Dex_{10kDa}$', r'$AHL_6$', r'$AHL_{12}$']
    xPos = [1,3,5,7]

    # print(diffu_df['diffCoef'])

    sns.boxplot(x='xPos', y="diffCoef", data=diffu_df, showfliers = False, hue = "mol", order=np.arange(8),
    palette='viridis', dodge=False, linewidth=0.5, ax = ax, width = 1)

    sns.stripplot(x='xPos', y="diffCoef", data=diffu_df,  hue = "mol", order=np.arange(8),
    palette='viridis', ax =ax, size = 2)

    ax.legend([],[], frameon=False)
    ax.grid(False)

    plt.xticks(xPos, xlabs, rotation = 30)
    ax.tick_params(axis='both', which='major', labelsize=6)

    # ax.set_xlabel('$Molecule$')
    ax.set_xlabel('')
    ax.set_ylabel(r'D $(\mu m^2/s)$')
    # ax.set_ylabel('Fold change (log2)')
    ax.set_xlim([0,8])
    ax.set_ylim(bottom=0)

    # trying to reduce the fontweight - doesn't work
    # ylabs = ax.get_yticklabels()
    # ax.set_xticklabels(xlabs, fontdict={'fontsize': 6,'fontweight': 100})
    # ax.set_yticklabels(ylabs, fontdict={'fontsize': 6,'fontweight': 100})

    plt.tight_layout()

    #remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    #
    figName_0 = f"diffuEsti.png"
    figName_1 = f"diffuEsti.pdf"
    # plt.show()
    plt.savefig(outPath/figName_0, transparent=True, dpi = 300)
    plt.savefig(outPath/figName_1, transparent=True)
    plt.close()

    print("analysis is complete.")
