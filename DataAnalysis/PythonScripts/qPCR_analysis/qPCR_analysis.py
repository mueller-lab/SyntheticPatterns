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
# from scipy.stats import linregress
import seaborn as sns
# import argparse

# functions

# main code
if __name__=="__main__":
    print("starting analysis ...")

    # global plot settings
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["font.size"] = 6
    plt.rcParams['legend.fontsize'] = 6
    plt.rcParams['pdf.fonttype'] = 42

    cwdPath=Path(os.path.abspath(os.getcwd()))
    inPath= cwdPath
    datetime_str=(datetime.now().strftime("%Y%m%d_"))
    # outDir = str(datetime_str+"output")
    outDir = "output"
    outPath=inPath/outDir
    outPath.mkdir(parents=True, exist_ok=True)

    # read data files and make one pandas df
    primers= ['LasR', 'mCherry', 'ihfB']
    strains = ['pAL201.5'] # +pAL201.5
    RT_NRT = ['RT','NRT'] # 1 = RT, 0 = NRT
    treats = [0]+[10**(i-1) for i in range(6) ] # AHL12 - 0, 10^-2 to 10^3 nM treatment for 4 h

    pri_ind = 7*(3*['LasR']+3*['mCherry']+3*['ihfB']+primers)+12*['none']
    rt_ind = 7*(9*['RT']+3*['NRT'])+12*['none']
    strn_ind= 96*['pAL201.5']

    treat_ind0 = [12*[treats[i]] for i in range(7)]
    treat_ind = [i for j in treat_ind0 for i in j]+12*['none']

    rep_ind = 8*(3*[1,2,3]+3*['none'])

    print(f'len {len(pri_ind), len(rt_ind), len(strn_ind), len(treat_ind), len(rep_ind)}')

    # Read data to DataFrame
    rawFile = inPath/"rawCq.xlsx"
    cq_df= pd.read_excel(str(rawFile),header=None)
    cq_df = cq_df.iloc[:,[5]] # remove excess rows and cols
    cq_df=cq_df.T
    cq_df.index = ['Cq']

    cq_df.columns =[pri_ind,rt_ind,strn_ind,treat_ind,rep_ind]
    cq_df.columns.set_names(['primer', 'RT', 'strain', 'treats', 'rep'], inplace=True)
    cq_df.to_csv(outPath/"Cq_data0.csv")

    cq_df = cq_df.drop(columns=['none'], level = 0)
    # cq_df = cq_df.drop(columns=['none'], level = 4)
    cq_df = cq_df.stack(level=[0,1,2,3,4])
    cq_df.index.set_names(['id', 'primer', 'RT', 'strain', 'treats', 'rep'], inplace=True)

    cq_df= cq_df.to_frame()
    cq_df.columns = ['cq_val']
    cq_df= cq_df.reset_index()
    cq_df = cq_df[['primer', 'RT', 'strain', 'treats', 'rep', 'cq_val']]
    cq_df = cq_df.set_index(['primer', 'RT', 'strain', 'treats', 'rep'])

    cq_df.to_csv(outPath/"Cq_data.csv")

    # process the data
    # steps - 1) calculate dCq_df 2) calculate ddCq_df 3) make boxplots

    for strn in strains:
        for pm in primers:
            gDf = cq_df.loc[(pm,'RT',strn)]
            refDF = cq_df.loc[('ihfB','RT',strn)]
            dCq_df = refDF - gDf

            dCq_df= dCq_df.loc[(treats),:]
            # print(dCq_df)

            dCq_df.columns=['dCq']
            dCq_csv = f'{strn}_{pm}_dCq.csv'
            dCq_df.to_csv(outPath/dCq_csv)

            dCq_C = dCq_df.loc[(0,[1,2,3]),'dCq'].mean()
            ddCq_df = dCq_df - dCq_C
            ddCq_df.columns=['ddCq']
            ddCq_df['foldCng']=2**(ddCq_df['ddCq'])
            ddCq_df['foldCng_log2']=ddCq_df['ddCq']

            # x-positions for the plot
            xScat = [1,1,1,3,3,3,5,5,5,7,7,7,9,9,9,11,11,11,13,13,13]
            ddCq_df['xPos']=xScat

            ddCq_csv = f'{strn}_{pm}_ddCq.csv'
            ddCq_df.to_csv(outPath/ddCq_csv)

            ddCq_df= ddCq_df.reset_index()

            #plot fold-change as box-plot with scatter
            fig, ax  = plt.subplots(nrows = 1, ncols = 1, figsize= (1.5,1.5))

            xlabs = [r'$0$', r'$10^{-2}$', r'$10^{-1}$', r'$10^{0}$', r'$10^{1}$', r'$10^{2}$', r'$10^{3}$']
            xPos = [1,3,5,7,9,11,13]

            sns.boxplot(x='xPos', y="foldCng", data=ddCq_df, showfliers = False, hue = "treats", order=np.arange(14),
            palette='magma', dodge=False, linewidth=0.5, ax = ax, width = 0.5)

            sns.stripplot(x='xPos', y="foldCng", data=ddCq_df,  hue = "treats", order=np.arange(14),
            palette='magma', ax =ax, size = 2)

            # sns.boxplot(x='xPos', y="foldCng_log2", data=ddCq_df, showfliers = False, hue = "treats", order=np.arange(14),
            # palette='magma', dodge=False, linewidth=1, ax = ax, width = 1)
            #
            # sns.stripplot(x='xPos', y="foldCng_log2", data=ddCq_df,  hue = "treats", order=np.arange(14),
            # palette='magma', ax =ax, size = 4)

            ax.legend([],[], frameon=False)
            ax.grid(False)

            plt.xticks(xPos, xlabs, rotation = 30)
            ax.tick_params(axis='both', which='major', labelsize=6)

            ax.set_xlabel(r'AHL$_{12}$ (nM)')
            ax.set_ylabel('Fold change')
            # ax.set_ylabel('Fold change (log2)')

            ax.set_xlim([0,14])

            yTop = np.max(ddCq_df['foldCng'])
            ax.set_ylim([0, max(8.0, yTop*1.1)])

            plt.tight_layout()
            #remove top and right spines
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            #
            figName_0 = f"{strn}_{pm}.png"
            figName_1 = f"{strn}_{pm}.pdf"
            # plt.show()
            plt.savefig(outPath/figName_0, transparent=True, dpi=300)
            plt.savefig(outPath/figName_1, transparent=True)
            plt.close()

    print("analysis is complete.")
