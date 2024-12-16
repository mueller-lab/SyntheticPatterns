"""
Amit Landge

Goal - To model positive loop circuit

"""

# import necessary packages
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
import time
import argparse
import os
from pathlib import Path
from datetime import datetime
import copy
from fipy import *
import pandas as pd
from scipy.optimize import curve_fit
import seaborn as sns


# define functions
class RD_simulate(object):
    """docstring for RD_simulate."""

    def __init__(self, nx, dx, dt, steps, vNames, resultPath, araList):
        super(RD_simulate, self).__init__()
        self.mesh= Grid1D(nx = nx, dx = dx)
        self.nx = nx
        self.dx = dx
        self.steps = steps
        self.dt = dt
        self.resultPath = resultPath
        self.vNames = vNames
        self.nodeN = len(vNames)
        self.doseRes = {}
        self.tAx = [i*dt for i in range(steps+2)]
        self.araList = araList

    def hillF(self,x,y,Kt, n = 2, k_leak = 0.01):
        """
        Hill function
        """
        return x*(k_leak + y**n/(Kt**n + y**n))

    def makeModel(self,y0,rates,D6):
        #unpack rates
        k_prot, kd_ssrA, k_sfGFP, kd_prot, k_ahl, kd_ahl, kd1_ahl, K_T, nH, k_leak, ara, K_Tara, nH_ara  = rates
        self.kd_ahl = kd_ahl
        self.kd1_ahl = kd1_ahl

        # make a separate directory for each arabinose conc results
        self.araPath=self.resultPath/f"ara_{ara:.2f}"
        self.araPath.mkdir(mode=0o777, parents=True, exist_ok=True)

        # save params to a text file
        with open(self.araPath/"params.csv", "w") as parFile:
            #write the header
            parFile.write("k_prot, kd_ssrA, k_sfGFP, kd_prot, k_ahl, kd_ahl, kd1_ahl, K_T, ara, nH, k_leak, K_Tara, nH_ara, D6")
            parFile.write('\n')
            parFile.write(f"{k_prot}, {kd_ssrA}, {k_sfGFP}, {kd_prot}, {k_ahl}, {kd_ahl}, {kd1_ahl}, {K_T}, {ara}, "+
            f"{nH}, {k_leak}, {K_Tara}, {nH_ara}, {D6}")
            parFile.write('\n')
        parFile.close()

        # starting state
        self.InVal= y0

        # create CellVariables
        self.vL = [] # list of CellVariables
        for i, vName in enumerate(self.vNames):
            var = CellVariable(name = vName, mesh = self.mesh, value = self.InVal[i])
            self.vL.append(var)

        #define the reaction part
        # note ---> ['LuxR', 'AHL6', 'sfGFP', 'AiiA', 'LuxI']
        self.LuxR_AHL6 = self.hillF(self.vL[0],self.vL[1],K_T,nH, k_leak)
        self.ara_indc = self.hillF(10.,ara,K_Tara,nH_ara,k_leak)

        self.eq00 =k_prot - kd_prot*self.vL[0] # LuxR
        self.eq01 =k_ahl*self.vL[4]  - self.vL[1]*(kd_ahl+kd1_ahl*self.vL[3]) # AHL6
        self.eq02 =k_sfGFP*self.LuxR_AHL6- kd_ssrA*self.vL[2] # sfGFP
        self.eq03 =k_prot*self.ara_indc- kd_ssrA*self.vL[3] # AiiA
        self.eq04 =k_prot*self.LuxR_AHL6- kd_ssrA*self.vL[4] # LuxI

        # define reaction-diffusion system of PDEs
        self.eq_LuxR = TransientTerm(var=self.vL[0]) == DiffusionTerm(coeff=0.0, var=self.vL[0])+ self.eq00
        self.eq_AHL6 = TransientTerm(var=self.vL[1]) == DiffusionTerm(coeff=D6, var=self.vL[1])+ self.eq01
        self.eq_sfGFP = TransientTerm(var=self.vL[2]) == DiffusionTerm(coeff=0.0, var=self.vL[2])+ self.eq02
        self.eq_AiiA = TransientTerm(var=self.vL[3]) == DiffusionTerm(coeff=0.0, var=self.vL[3])+ self.eq03
        self.eq_LuxI = TransientTerm(var=self.vL[4]) == DiffusionTerm(coeff=0.0, var=self.vL[4])+ self.eq04

        self.eq = self.eq_LuxR & self.eq_AHL6 & self.eq_sfGFP & self.eq_AiiA & self.eq_LuxI

        #write results to .csv files
        for i, vName in enumerate(self.vNames):
            fName = f"{vName}.csv"
            with open(self.araPath/fName, "w") as parFile:
                #write the header
                parFile.write("t_min")
                for j in range(self.nx):
                    parFile.write(f",{j*self.dx}")
                parFile.write('\n')

                #write the first row
                parFile.write("0.00")
                for j in range(self.nx):
                    parFile.write(f",{self.vL[i][j]}")
                parFile.write('\n')
            parFile.close()

    def simulateModel(self, ):
        # iterate to numerically solve the system and save result to a dictionary
        step = 0
        while step <= self.steps:
            self.eq.solve(dt=self.dt)

            #save results to the .csv files
            for j, vName in enumerate(self.vNames):
                fName = f"{vName}.csv"
                with open(self.araPath/fName, "a") as parFile:
                    #append a new row to the file
                    parFile.write(f"{self.dt*(step+1):.2f}")
                    for k in range(self.nx):
                        parFile.write(f",{self.vL[j][k]}")
                    parFile.write('\n')
                parFile.close()
            step+=1

    def analysis_plots(self, subVal = 4, xticInt = 60):
        """
        method to -
        1. make space-time plots (not anymore)
        2. save dose-response data to a dict
        """
        # dict to store time series of all variables
        timeSerDict={}

        # plot the results as space-time graphs and save
        for j, vName in enumerate(self.vNames):
            fName = f"{vName}.csv"
            fPath = self.araPath/fName
            data_df = pd.read_csv(str(fPath),header=0, index_col = 0)

            vals = data_df.to_numpy(dtype='float')
            vals_sub = vals[::subVal,:]

            # get an array of mean values at the center (10% of the total region)
            meanCent = np.mean(vals[:,int(0.4*self.nx):int(0.5*self.nx)], axis = 1)

            # add timeseries to the dict
            timeSerDict[vName]= meanCent

            # add final sfGFP conc. to a dict to make the dose response curve
            if vName =="sfGFP" and self.tAx[-1]>=360:
                self.doseRes[ara] =meanCent[int(360/self.dt)]

            if vName =="sfGFP" and self.tAx[-1]<360:
                self.doseRes[ara] =meanCent[-1]

        # save time series of all variable in the respective ara dir
        timeSerDF = pd.DataFrame.from_dict(timeSerDict)
        timeSerDF.to_csv(self.araPath/'timeSerAll.csv')

    def timeSerPlots(self, ):
        """
        to make time-series plots of the variable output
        """
        # make a color dict to assign colors to typeId
        viridis = matplotlib.colormaps.get_cmap('viridis', 256)
        cmap0 = viridis(np.linspace(0, 1, len(self.araList)))

        for vName in self.vNames:
            # load sfGFP timeseries at all ahl treatments
            normSeries = []
            for ara in self.araList:
                fName = self.resultPath/f"ara_{ara:.2f}"/'timeSerAll.csv'
                data_df = pd.read_csv(str(fName),header=0, index_col = 0)
                vals = np.array(data_df[vName])
                normSeries.append(vals)

            normArr = np.array(normSeries)

            # make plots without scaling the data to maximum of 1
            fig, ax  = plt.subplots(nrows = 1, ncols = 1, figsize= (2.4,1.5))

            for i, ara in enumerate(self.araList):
                ax.plot(self.tAx, normArr[i,:], linewidth = 1, c=cmap0[i], label = f'{ara:.0e} nM')

            ax.set_xlabel("Time (min)")
            ax.set_ylabel("Simulated sfGFP output")
            ax.set_ylim(0.0,1.02*np.amax(normArr))
            ax.set_xlim(0.0,1.02*self.tAx[-1])
            ax.tick_params(axis='both', which='major', labelsize=6)
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title = "Arabinose")

            #remove top and right spines
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            plt.tight_layout()

            figName_0 = f"TimeSeries_{vName}.png"
            figName_1 = f"TimeSeries_{vName}.pdf"

            # plt.show()
            plt.savefig(self.resultPath/figName_0, transparent=True, dpi = 300)
            plt.savefig(self.resultPath/figName_1, transparent=True)
            plt.close()

            # make plots after scaling the data to maximum of 1
            normArr0 = normArr/np.amax(normArr)

            # to make time-series plots with mean values from domain center
            fig, ax  = plt.subplots(nrows = 1, ncols = 1, figsize= (2.4,1.5))

            for i, ara in enumerate(self.araList):
                ax.plot(self.tAx, normArr0[i,:], linewidth = 1, c=cmap0[i], label = f'{ara:.0e} nM')

            ax.set_xlabel("Time (min)")
            ax.set_ylabel("Simulated sfGFP output")
            ax.set_ylim(0.0,1.02)
            ax.set_xlim(0.0,1.02*self.tAx[-1])
            ax.tick_params(axis='both', which='major', labelsize=6)
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title = "Arabinose")

            #remove top and right spines
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            plt.tight_layout()

            figName_0 = f"normTimeSeries_{vName}.png"
            figName_1 = f"normTimeSeries_{vName}.pdf"

            # plt.show()
            plt.savefig(self.resultPath/figName_0, transparent=True, dpi = 300)
            plt.savefig(self.resultPath/figName_1, transparent=True)
            plt.close()

    def makeDoseResponse(self,):
        """
        method to make dose-response plots with mean values from domain center
        """
        fig, ax  = plt.subplots(nrows = 1, ncols = 1, figsize= (1.5,1.5))

        # make data array
        doseRes = np.array(list(self.doseRes.items()))

        #save doseRes dict as a pandas DataFrame
        doseDf = pd.DataFrame(self.doseRes.items(), columns = ['Ara', 'sfGFP'])
        doseDf.to_csv(self.resultPath/"doseRes.csv")

        xAra = doseRes[1:,0]
        yRes = doseRes[1:,1]/np.amax(doseRes[1:,1])
        ax.scatter(xAra, yRes, marker= "o", c='k', s = 5)

        ####
        # dict to store fit parameters
        fitDict = {"K":[], "n":[], "err_K":[], "err_n":[], "ymin":[], "ymax":[]}

        #fit data to rev_sigmoidFunc - 1) use the x and yAvg to calculate the fit parameters; 2) define to new xFit array and calculate yFit using the
        # obtained fit parameters- opt_K and opt_n
        bnds= ([10**1, 0.5],[10**5, 2.0])
        mthd = 'trf'

        xdata=np.zeros((len(xAra), 3))
        xdata[:,0]= xAra
        xdata[:,1]= np.amin(yRes)
        xdata[:,2]= np.amax(yRes)

        popt, pcov = curve_fit(self.rev_sigmoidFunc, xdata, yRes, bounds=bnds,  method=mthd)
        perr = np.sqrt(np.diag(pcov))

        xFit = np.logspace(0,6,50)
        opt_K, opt_n = popt

        fitDict["K"].append(opt_K)
        fitDict["n"].append(opt_n)
        fitDict["err_K"].append(perr[0])
        fitDict["err_n"].append(perr[1])
        fitDict["ymin"].append(np.amin(yRes))
        fitDict["ymax"].append(np.amax(yRes))

        xdata=np.zeros((len(xFit), 3))
        xdata[:,0]= xFit
        xdata[:,1]= np.amin(yRes)
        xdata[:,2]= np.amax(yRes)

        yFit = self.rev_sigmoidFunc(xdata,opt_K,opt_n)

        ax.plot(xFit,yFit, c="k", linewidth = 1)
        ####
        # ax.plot(xAra, yRes, 'g--', linewidth = 1)

        ax.set_xscale("log")
        ax.set_xlabel("Arabinose (nM)")
        ax.set_ylabel("Simulated sfGFP output")
        ax.set_ylim(0.0,1.02*np.amax(yRes))
        ax.set_xlim(xAra[0],1.2*xAra[-1])

        #remove top and right spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(False)
        plt.tight_layout()

        figName_0 = f"ara_doseRe.png"
        figName_1 = f"ara_doseRe.pdf"

        # plt.show()
        plt.savefig(self.resultPath/figName_0, transparent=True, dpi = 300)
        plt.savefig(self.resultPath/figName_1, transparent=True)
        plt.close()

        fit_df = pd.DataFrame(fitDict,index=[360])
        fit_df.to_csv(self.resultPath/f"fitParams.csv")

    def rev_sigmoidFunc(self, xdata, K, n):
        x = xdata[:,0]
        ymin=xdata[:,1]
        ymax = xdata[:,2]
        y = ymin + (ymax-ymin)*K**n/(x**n+K**n)
        return y

# main code
if __name__=="__main__":
    # get current dir and create result dir
    print("\n**************************************************** \nStarting the program...")
    t1 = time.time()

    # global plot settings
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["font.size"] = 6
    plt.rcParams['legend.fontsize'] = 6
    plt.rcParams['pdf.fonttype'] = 42

    #get current dir and datetime
    cwdPath=Path(os.path.abspath(os.getcwd()))
    now=datetime.now()
    datetime_str=(now.strftime("%Y%m%d_")) # %H%M%S_

    #make output directory
    dirName=str(datetime_str+ "mod3_posLoop") # the positive loop model

    resultPath=cwdPath/'data_1D'/dirName
    resultPath.mkdir(mode=0o777, parents=True, exist_ok=True)
    print(f"Created result directory {dirName} at {time.time() - t1} sec ...")

    #define global parameters
    steps = 720
    dt=1.0 # time in min
    gridSize = 2 # mm
    dx= 0.1 # mm
    nx = int(gridSize/dx)
    vNames = ['LuxR', 'AHL6', 'sfGFP', 'AiiA', 'LuxI']

    # define diffusion rates (per min)
    D6= 0. # 3.6*1E-2 # AHL6 diffusion rate(act)

    # synthesis rates (1/min, )
    k_prot = 0.2 # protein synthesis (transcrition + translation)
    k_sfGFP = 0.1 # sfGFP protein synthesis + maturation
    k_ahl = 0.002 # AHL synthesis rate (by LuxI)

    # degradation rates (unit = 1/min, )
    kd_ssrA = 0.018 # sfGFP decay (LVA tagged)
    kd_prot = 0.0028 # LuxR decay, half-life of 4 h
    kd_ahl = 0.0011 # ahl decay (hydroslysis), half-life ~ 10 h
    kd1_ahl = 0.01 # ahl degradation by lactonase (arabinose-induced)

    # activation threshold
    K_T = 10.0 # nM ahl6 threshold conc.
    K_Tara = 6000. # arabinose induction threshold
    nH = 1.5 # Hill coefficient
    nH_ara = 1.5
    k_leak = 0.08 # leaky expression of sfGFP from pLux promoter (same leakiness for pBAD)

    araList = [0., 1.,10.,100., 1000.,1E4, 1E5, 1E6] # 1.,10.,100.,
    # araList = [0., 1., 1000., 1E5] # 1.,10.,100.,

    mod_posLoop = RD_simulate(nx, dx, dt, steps, vNames, resultPath, araList)

    # y0 = [1.0 for i in range(3)]
    y0 = [60.,1.0, 180.0,  0.10, 360.0] # 'LuxR', 'AHL6', 'sfGFP', 'AiiA', 'LuxI',
    # y0 = 5*[40.] # 'LuxR', 'AHL6', 'sfGFP', 'AiiA', 'LuxI',

    for ara in araList:
        rates = [k_prot, kd_ssrA, k_sfGFP, kd_prot, k_ahl, kd_ahl, kd1_ahl, K_T, \
        nH, k_leak, ara, K_Tara, nH_ara]
        mod_posLoop.makeModel(y0,rates,D6)
        mod_posLoop.simulateModel()
        mod_posLoop.analysis_plots()
        print(f'Completed ara dose = {ara}')

    mod_posLoop.makeDoseResponse()
    mod_posLoop.timeSerPlots()

    #end code
    t2= time.time()
    totalSec= t2-t1
    Sec=int(totalSec%60)
    Hrs=int(totalSec//3600)
    Min=int((totalSec%3600)//60)

    print ("Program completed in %sHr:%sMin:%ssec\n"%(Hrs,Min,Sec))
