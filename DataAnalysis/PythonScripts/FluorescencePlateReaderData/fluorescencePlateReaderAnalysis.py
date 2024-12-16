"""
Amit Landge
Custom code to analyse and plot data from the plate reader assays.

"""
import os
import time
from datetime import datetime
from pathlib import Path
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from scipy.optimize import curve_fit
import seaborn as sns
# import argparse

# functions

def rev_sigmoidFunc(xdata, K, n):
    x = xdata[:,0]
    ymin=xdata[:,1]
    ymax = xdata[:,2]
    y = ymin + (ymax-ymin)*K**n/(x**n+K**n)
    return y

# main code
if __name__=="__main__":
    print("starting analysis ...")

    # global plot settings
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["font.size"] = 6
    plt.rcParams['legend.fontsize'] = 6
    plt.rcParams['pdf.fonttype'] = 42

    cwdPath=Path(os.path.abspath(os.getcwd()))
    #get paths to read data from each biological replicate
    inPaths= [cwdPath/'20230628_posLoopAll', \
    cwdPath/'20230629_posLoopAll', \
    cwdPath/'20230630_posLoopAll']

    datetime_str=(datetime.now().strftime("%Y%m%d_")) # %H%M%S_
    # outDir = str(datetime_str+"outputPlateReaderAnalysis")
    outDir = "OutputPlateReaderAnalysis"
    outPath=cwdPath/outDir
    outPath.mkdir(parents=True, exist_ok=True)

    # read data files and make one pandas df
    timePoints = [5*i for i in range(145)] # in minutes
    readChanls = ["OD600", "sfGFP"]
    strains = 8*(2*["NoPlasmid"]+2*["pAL101"]+2*["pAL102"]+2*["pAL103"]+2*["pAL104"]+2*["NoCells"])

    treatments0= [0,1, 10**1, 10**2, 10**3, 10**4, 10**5, 10**6]
    treatments1 = [12*[trt] for trt in treatments0]
    treatments = [i for j in treatments1 for i in j]

    tech_reps = 8*(6*[1,2])


    # for each biological replicate, read data to DataFrame

    # perform background subtraction while reading each biological replicate data
    OD_back = 3*[0.09]
    sfGFP_back = [0, 0, 0]

    OD_DataDfs = {} # to hold the DFs
    sfGFP_DataDfs = {} #
    for i, inPath in enumerate(inPaths):
        # Read data to DataFrame
        sfGFP_File = inPath/"sfGFP.xlsx"
        OD_File= inPath/"od600.xlsx"
        sfGFP_df= pd.read_excel(str(sfGFP_File),header=0)
        OD_df= pd.read_excel(str(OD_File),header=0)

        sfGFP_df= sfGFP_df.T
        OD_df= OD_df.T

        sfGFP_df = sfGFP_df.iloc[1:146,2:] # remove excess rows and cols
        OD_df = OD_df.iloc[1:146,2:] # remove excess rows and cols

        #index as per the time in minutes
        sfGFP_df.index= timePoints
        OD_df.index= timePoints

        # name the columns
        sfGFP_df.columns= [strains, treatments, tech_reps]
        sfGFP_df.columns.set_names(['strain', 'treat', 'rep'], inplace=True)
        # sfGFP_df.drop(columns=['none'], inplace=True)
        sfGFP_df= sfGFP_df.T

        # sfGFP_df = sfGFP_df.stack(level=[0,1,2])
        # sfGFP_df.index.set_names(['strain', 'treat', 'rep'], inplace=True)
        # print(sfGFP_df.head)

        OD_df.columns= [strains, treatments, tech_reps]
        OD_df.columns.set_names(['strain', 'treat', 'rep'], inplace=True)
        # OD_df.drop(columns=['none'], inplace=True)
        OD_df= OD_df.T

        #background subtraction
        sfGFP_df = sfGFP_df - sfGFP_back[i]
        OD_df = OD_df - OD_back[i]

        # add df to dict
        sfGFP_DataDfs[i+1] = sfGFP_df
        OD_DataDfs[i+1] = OD_df

    # concatenate the biological replicate DFs to make one big DF
    sfGFP_df = pd.concat(sfGFP_DataDfs, axis = 0, names = ['Replicate'])
    OD_df = pd.concat(OD_DataDfs, axis = 0, names = ['Replicate'])

    #save raw data to .csv files
    OD_df.to_csv(outPath/"OD_data.csv")
    sfGFP_df.to_csv(outPath/"sfGFP_data.csv")

    print("cleaning and labeling is complete.")

    # process the data - 1) Background subtraction-->  sfGFP - sfGFP_back,and OD - OD_back;
    # 2) Normalize sfGFP signal by OD; 3) Make plots

    """
    Idea -
    The autofluorescence signal changes over time, this is seen in the NoPlasmid controls as well as
    the enigineered E. coli strains are also affected by it.
    To correct for this - The NoPlasmid autofluorescence signal for sfGFP and mCherry should be
    subtracted from the corresponding signal for the other strains.
    The new data will be stored as sfGFP_minusAuto.csv and mCherry_minusAuto.csv
    """

    # list the strains
    strns = ['NoPlasmid', 'pAL101', 'pAL102', 'pAL103', 'pAL104'] # 'NoCells'

    sfGFP_df = sfGFP_df.reset_index()

    NoPlasmid_sfGFP = sfGFP_df.loc[sfGFP_df['strain']=='NoPlasmid']
    NoPlasmid_sfGFP.index = [i for i in range(NoPlasmid_sfGFP.shape[0]) ]

    sfGFP_mAU_Dict = {}
    for strn in strns:
        sfGFP_mAU = sfGFP_df.loc[sfGFP_df['strain']==strn][timePoints]
        sfGFP_mAU.index = [i for i in range(sfGFP_mAU.shape[0])]

        sfGFP_mAU = sfGFP_mAU - NoPlasmid_sfGFP[timePoints]
        toJoin_sfGFP_df = sfGFP_df.loc[sfGFP_df['strain']==strn][['Replicate', 'treat', 'rep']]
        toJoin_sfGFP_df.index =  [i for i in range(sfGFP_mAU.shape[0])]

        sfGFP_mAU = toJoin_sfGFP_df.join(sfGFP_mAU)

        sfGFP_mAU_Dict[strn] = sfGFP_mAU

    # print(sfGFP_mAU_Dict['pAL101'].head)

    # concatenate corrected data by strains
    sfGFP_mAU_df = pd.concat(sfGFP_mAU_Dict, axis = 0, names = ['strain'])
    sfGFP_mAU_df = sfGFP_mAU_df.reset_index()
    sfGFP_mAU_df = sfGFP_mAU_df.set_index(['Replicate','treat', 'strain', 'rep'])
    sfGFP_mAU_df = sfGFP_mAU_df.drop(columns = 'level_1')

    #save to .csv files
    sfGFP_mAU_df.to_csv(outPath/"sfGFP_mAU_data.csv")

    # Background subtraction was performed above
    sfGFP_df_backSub = sfGFP_mAU_df
    OD_df_backSub = OD_df

    # print(OD_df_backSub.head)

    sfGFP_by_OD_df = sfGFP_df_backSub / OD_df_backSub
    sfGFP_by_OD_df.reset_index(inplace=True) # to turn index into selectable columns
    sfGFP_by_OD_df.to_csv(outPath/"sfGFP_by_OD_data.csv")

    OD_df_backSub.reset_index(inplace=True)
    sfGFP_df_backSub.reset_index(inplace=True)

    # get data for the pBC strain
    strns = ['pAL101', 'pAL102', 'pAL103', 'pAL104'] # 'NoCells'
    for strn in strns:
        data_pBC = sfGFP_by_OD_df.loc[sfGFP_by_OD_df['strain'] == strn]

        # only select the Ara treated samples (drop rows - "Control" and "none")
        # data_pBC = data_pBC.set_index(["treat"])
        #
        # data_pBC = data_pBC.drop(index = ["Control"])
        # data_pBC = data_pBC.sort_index()
        #
        # data_pBC = data_pBC.reset_index()

        data_pBC = data_pBC.set_index(['Replicate','treat', 'strain', 'rep'])
        # print(data_pBC.head)

        #normalize the data to scale between 0 to 1
        allMin = np.amin(data_pBC.min())
        allMax = np.amax(data_pBC.max())
        # print(f"for pBC; allMin={allMin}, allMax={allMax}")

        data_norm = data_pBC / allMax
        data_norm =data_norm.reset_index()
        data_norm.to_csv(outPath/f"{strn}_norm_data.csv")

        # make dose-response plots at every time-point
        # dict to store fit parameters
        fitDict = {"K":[], "n":[], "err_K":[], "err_n":[], "ymin":[], "ymax":[]}

        doseTPs = [120, 360] # in minutes
        for t0 in doseTPs: # from 4 h to 8 h --> 240 to 480 min -> 48 to 96
            #select the data
            data01 = data_norm[['Replicate','treat','rep', t0]]

            # drop 'treat' = 0
            data01 = data01.loc[data01['treat']!=0]

            # add a column with normalized values
            normT0 = data01[t0]/np.max(data01[t0])
            data01.insert(1, 'normT0', normT0)

            # print(f'data01 is \n {data01}')

            # make array of the Ara concentrations
            x = np.unique(np.array(data01['treat']))

            # make dose-response plot
            fig, ax  = plt.subplots(nrows = 1, ncols = 1, figsize= (1.5,1.5))

            #make errorbar (standard deviation)
            # ydata = ydata.astype(float)
            # err = np.std(ydata, axis = 0)
            # yAvg = np.mean(ydata, axis = 0)
            err= data01.groupby('treat', sort= False).agg({'normT0':'std'})['normT0']
            yAvg= data01.groupby('treat', sort= False).agg({'normT0':'mean'})['normT0']

            err= np.array(err)
            yAvg= np.array(yAvg)

            # print(f'yAvg is\n {yAvg}')

            #fit data to rev_sigmoidFunc - 1) use the x and yAvg to calculate the fit parameters; 2) define to new xFit array and calculate yFit using the
            # obtained fit parameters- opt_K and opt_n
            bnds= ([10**2, 1.0],[10**5, 3.0])
            mthd = 'trf'

            xdata=np.zeros((len(x), 3))
            xdata[:,0]= x
            xdata[:,1]= np.amin(yAvg)
            xdata[:,2]= np.amax(yAvg)

            popt, pcov = curve_fit(rev_sigmoidFunc, xdata, yAvg, bounds=bnds,  method=mthd)
            perr = np.sqrt(np.diag(pcov))

            xFit = np.logspace(0,6,50)
            opt_K, opt_n = popt

            fitDict["K"].append(opt_K)
            fitDict["n"].append(opt_n)
            fitDict["err_K"].append(perr[0])
            fitDict["err_n"].append(perr[1])
            fitDict["ymin"].append(np.amin(yAvg))
            fitDict["ymax"].append(np.amax(yAvg))

            xdata=np.zeros((len(xFit), 3))
            xdata[:,0]= xFit
            xdata[:,1]= np.amin(yAvg)
            xdata[:,2]= np.amax(yAvg)

            yFit = rev_sigmoidFunc(xdata,opt_K,opt_n)

            ax.plot(xFit,yFit, c="k", linewidth = 1)
            markers, caps, bars = ax.errorbar(x,yAvg, yerr= err, c ='k', fmt="none", capsize = 1.0, capthick = 1.0)
            [bar.set_alpha(0.2) for bar in bars] # to make the errorbars transparent
            [cap.set_alpha(0.2) for cap in caps]

            # for i, yRep in enumerate(ydata):
                # ax.scatter(x, yRep, marker= ".", c='k', s = 20)
            # add scatter
            # ax.scatter(data01['treat'], data01['normT0'], marker= ".", c='g', s = 10)
            # sns.scatterplot(data=data01, x="treat", y="normT0", hue="Replicate", edgecolor="none", s = 10)
            sns.scatterplot(data=data01, x="treat", y="normT0", style="Replicate", \
            edgecolor="none", s = 10, hue="Replicate", palette = ['k', 'k', 'k']) # markers = ['o','X',],

            ax.set_xscale("log")
            ax.grid(False)

            ax.set_xlabel("Arabinose (nM)")
            ax.set_ylabel(r"Normalized fluorescence")
            # ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

            ax.set_ylim(0,1.02)
            ax.set_xlim(x[0], x[-1]*1.2)

            # ax.text(x[-3], yAvg[1], f" K= {opt_K:.1f} nM")

            ax.tick_params(axis='both', which='major', labelsize=6)

            plt.tight_layout()

            #remove top and right spines
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

            figName_0 = f"{strn}_time_{t0}_doseRe.png"
            figName_1 = f"{strn}_time_{t0}_doseRe.pdf"
            # plt.show()
            plt.savefig(outPath/figName_0, transparent=True, dpi = 300.)
            plt.savefig(outPath/figName_1, transparent=True)
            plt.close()

        fit_df = pd.DataFrame(fitDict,index=doseTPs)
        fit_df.to_csv(outPath/f"{strn}_fitParams.csv")

###########################################################################################
        # make timeseries plots for all AHL treatments (all treats in one plot color coded)
        fig, ax  = plt.subplots(nrows = 1, ncols = 1, figsize= (2.4,1.5))
        x = timePoints

        treatments = np.unique(np.array(data_norm.loc[data_norm['rep']==1]['treat']))
        treatments = np.sort(treatments)

        # make a color dict to assign colors to typeId
        viridis = matplotlib.colormaps.get_cmap('viridis')
        cmap0 = viridis(np.linspace(0, 1, len(treatments)))

        # plot a timeseries of fluorescence/OD_600
        for j, treat in enumerate(treatments):
            treatColor = cmap0[j]
            y = np.array(data_norm.loc[data_norm['treat']==treat][timePoints])
            # print(y)

            #make errorbar (standard deviation)
            y = y.astype(float)
            err = np.std(y, axis = 0)
            yAvg = np.mean(y, axis = 0)

            # markers, caps, bars = ax.errorbar(x,yAvg, yerr= err, fmt="none", capsize = 0.5, capthick= 1.0, c = treatColor)
            # [bar.set_alpha(0.2) for bar in bars] # to make the errorbars transparent
            # [cap.set_alpha(0.2) for cap in caps]

            plt.fill_between(x,yAvg-err,yAvg+err, alpha = 0.2, color = treatColor, edgecolor= 'face')

            ax.plot(x, yAvg, linewidth=1, c=treatColor, label = f"{treat:.0e} nM")

            # for i, yRep in enumerate(y):
            #     ax.scatter(x, yRep, marker= ".", color=treatColor, s = 20)

        ax.set_xlabel("Time (min)")
        ax.set_ylabel(r"Normalized fluorescence")
        ax.set_ylim(0,1.02)
        ax.set_xlim(0, x[-1]*1.02)
        ax.tick_params(axis='both', which='major', labelsize=6)
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

        #remove top and right spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        plt.tight_layout()
        figName_0 = f"{strn}_timeSeries.png"
        figName_1 = f"{strn}_timeSeries.pdf"
        # plt.show()
        plt.savefig(outPath/figName_0, transparent=True, dpi = 300.)
        plt.savefig(outPath/figName_1, transparent=True)
        plt.close()

#######################################
        # plot a timeseries of OD_600
        data_norm = OD_df_backSub.loc[OD_df_backSub['strain'] == strn]

        # only select the Ara treated samples (drop rows - "Control" and "none")
        # data_norm = data_norm.set_index(["treat"])
        #
        # data_norm = data_norm.drop(index = ["Control"])
        # data_norm = data_norm.sort_index()
        #
        # data_norm = data_norm.reset_index()

        fig, ax  = plt.subplots(nrows = 1, ncols = 1, figsize= (2.4,1.5))
        x = timePoints

        treatments = np.unique(np.array(data_norm.loc[data_norm['rep']==1]['treat']))
        treatments = np.sort(treatments)

        # make a color dict to assign colors to typeId
        viridis = matplotlib.colormaps.get_cmap('viridis')
        cmap0 = viridis(np.linspace(0, 1, len(treatments)))

        for j, treat in enumerate(treatments):
            treatColor = cmap0[j]
            y = np.array(data_norm.loc[data_norm['treat']==treat][timePoints])
            # print(y)

            #make errorbar (standard deviation)
            y = y.astype(float)
            err = np.std(y, axis = 0)
            yAvg = np.mean(y, axis = 0)

            # markers, caps, bars = ax.errorbar(x,yAvg, yerr= err, fmt="none", capsize = 0.5, capthick= 1.0, c = treatColor)
            # [bar.set_alpha(0.2) for bar in bars] # to make the errorbars transparent
            # [cap.set_alpha(0.2) for cap in caps]

            plt.fill_between(x,yAvg-err,yAvg+err, alpha = 0.2, color = treatColor, edgecolor= 'face')

            ax.plot(x, yAvg, linewidth=1, c=treatColor, label = f"{treat:.0e} nM")

            # for i, yRep in enumerate(y):
            #     ax.scatter(x, yRep, marker= ".", color=treatColor, s = 20)

        ax.set_xlabel("Time (min)")
        ax.set_ylabel(r"OD$_{600}$ (A.U.)")
        ax.set_ylim(0,)
        ax.set_xlim(0, x[-1]*1.02)
        ax.tick_params(axis='both', which='major', labelsize=6)
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

        #remove top and right spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        plt.tight_layout()
        figName_0 = f"{strn}_timeSeries_OD.png"
        figName_1 = f"{strn}_timeSeries_OD.pdf"
        # plt.show()
        plt.savefig(outPath/figName_0, transparent=True, dpi = 300.)
        plt.savefig(outPath/figName_1, transparent=True)
        plt.close()

##############################################
        # plot a timeseries of fluorescence
        data_norm = sfGFP_df_backSub.loc[sfGFP_df_backSub['strain'] == strn]

        # only select the Ara treated samples (drop rows - "Control" and "none")
        # data_norm = data_norm.set_index(["treat"])
        #
        # data_norm = data_norm.drop(index = ["Control"])
        # data_norm = data_norm.sort_index()
        #
        # data_norm = data_norm.reset_index()

        fig, ax  = plt.subplots(nrows = 1, ncols = 1, figsize= (2.4,1.5))
        x = timePoints

        treatments = np.unique(np.array(data_norm.loc[data_norm['rep']==1]['treat']))
        treatments = np.sort(treatments)

        # make a color dict to assign colors to typeId
        viridis = matplotlib.colormaps.get_cmap('viridis')
        cmap0 = viridis(np.linspace(0, 1, len(treatments)))

        for j, treat in enumerate(treatments):
            treatColor = cmap0[j]
            y = np.array(data_norm.loc[data_norm['treat']==treat][timePoints])
            # print(y)

            #make errorbar (standard deviation)
            y = y.astype(float)
            err = np.std(y, axis = 0)
            yAvg = np.mean(y, axis = 0)
            # markers, caps, bars = ax.errorbar(x,yAvg, yerr= err, fmt="none", capsize = 0.5, capthick= 1.0, c = treatColor)
            # [bar.set_alpha(0.2) for bar in bars] # to make the errorbars transparent
            # [cap.set_alpha(0.2) for cap in caps]

            plt.fill_between(x,yAvg-err,yAvg+err, alpha = 0.2, color = treatColor, edgecolor= 'face')

            ax.plot(x, yAvg, linewidth=1, c=treatColor, label = f"{treat:.0e} nM")

            # for i, yRep in enumerate(y):
            #     ax.scatter(x, yRep, marker= ".", color=treatColor, s = 20)

        ax.set_xlabel("Time (min)")
        ax.set_ylabel("Fluorescence (A.U.)")
        ax.set_ylim(0,)
        ax.set_xlim(0, x[-1]*1.02)
        ax.tick_params(axis='both', which='major', labelsize=6)
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

        #remove top and right spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        plt.tight_layout()
        figName_0 = f"{strn}_timeSeries_Fluor.png"
        figName_1 = f"{strn}_timeSeries_Fluor.pdf"
        # plt.show()
        plt.savefig(outPath/figName_0, transparent=True, dpi = 300.)
        plt.savefig(outPath/figName_1, transparent=True)
        plt.close()

    print("analysis is complete.")
