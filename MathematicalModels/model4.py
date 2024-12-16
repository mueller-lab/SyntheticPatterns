"""
Amit Landge

Goal - To model negative loop circuit
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

    def __init__(self, nx, dx, dt, steps, vNames, resultPath, ahl12List):
        super(RD_simulate, self).__init__()
        self.mesh= Grid1D(nx = nx, dx = dx)
        self.nx = nx
        self.dx = dx
        self.tAx = [i*dt for i in range(steps+2)]
        self.steps = steps
        self.dt = dt
        self.resultPath = resultPath
        self.vNames = vNames
        self.nodeN = len(vNames)
        self.doseRes = {}
        self.ahl12List = ahl12List

    def hillF(self,x,y,Kt, n = 2, k_leak = 0.01):
        """
        Hill function
        """
        return x*(k_leak + y**n/(Kt**n + y**n))

    def makeModel(self,y0,rates,D6,D12):
        """
        define the system of equations
        """
        #unpack rates
        k_prot, kd_prot, kd_ssrA, k_mChr, K_T6, K_T12, K_T12xR, K_Taqs, nH_6, nH_12, nH_aqs, LuxR, ahl12, AHL6, k_leak  = rates

        # make a separate directory for each ahl12 conc results
        self.ahl12Path=self.resultPath/f"ahl12_{ahl12:.3f}"
        self.ahl12Path.mkdir(mode=0o777, parents=True, exist_ok=True)

        # save params to a text file
        with open(self.ahl12Path/"params.csv", "w") as parFile:
            #write the header
            parFile.write("k_prot, kd_prot, kd_ssrA, k_mChr, K_T6, \
            K_T12, K_T12xR, K_Taqs, nH_6, nH_12, nH_aqs, LuxR, ahl12, AHL6, k_leak")
            parFile.write('\n')
            parFile.write(f"{k_prot}, {kd_prot}, {kd_ssrA}, {k_mChr},{K_T6},\
            {K_T12},{K_T12xR},{K_Taqs}, {nH_6}, {nH_12}, {nH_aqs}, {LuxR}, {ahl12},{AHL6},{k_leak}")
            parFile.write('\n')
        parFile.close()

        # starting state
        self.InVal= y0

        # create CellVariables
        self.vL = [] # list of CellVariables
        for i, vName in enumerate(self.vNames):
            var = CellVariable(name = vName, mesh = self.mesh, value = self.InVal[i])
            self.vL.append(var)

        #define the reaction part - ['LasR', 'mCherry', 'Aqs1']
        self.LuxR_AHL6 = self.hillF(LuxR, AHL6, K_T6, nH_6, k_leak)
        self.LuxR_AHL12 = self.hillF(LuxR, ahl12, K_T12xR, nH_12, k_leak)

        self.LasR_inh = self.hillF(1.0, K_Taqs, self.vL[2], nH_aqs,k_leak)
        self.LasR_AHL12 = self.hillF(self.vL[0], ahl12, K_T12/self.LasR_inh, nH_12, 0.02)

        self.eq00 = k_prot - kd_prot*self.vL[0] # LasR
        self.eq01 = k_mChr*self.LasR_AHL12 - kd_ssrA*self.vL[1] # mCherry
        self.eq02 = k_prot*(self.LuxR_AHL6+self.LuxR_AHL12)- kd_prot*self.vL[2] # Aqs1

        self.eq_LasR = TransientTerm(var=self.vL[0]) == DiffusionTerm(coeff=0.0, var=self.vL[0])+ self.eq00
        self.eq_mChr = TransientTerm(var=self.vL[1]) == DiffusionTerm(coeff=0.0, var=self.vL[1])+ self.eq01
        self.eq_aqs = TransientTerm(var=self.vL[2]) == DiffusionTerm(coeff=0.0, var=self.vL[2])+ self.eq02

        self.eq = self.eq_LasR & self.eq_mChr & self.eq_aqs

        self.mid0 = int(0.4*self.nx)
        self.mid1 = int(0.6*self.nx)

        #write results to .csv files
        for i, vName in enumerate(self.vNames):
            fName = f"{vName}.csv"
            with open(self.ahl12Path/fName, "w") as parFile:
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
        """
        use eq.solve to interatively solve the system
        """
        # iterate to numerically solve the system and save result to a dictionary
        step = 0
        while step <= self.steps:
            self.eq.solve(dt=self.dt)

            #save results to the .csv files
            for j, vName in enumerate(self.vNames):
                fName = f"{vName}.csv"
                with open(self.ahl12Path/fName, "a") as parFile:
                    #append a new row to the file
                    parFile.write(f"{self.dt*(step+1):.2f}")
                    for k in range(self.nx):
                        parFile.write(f",{self.vL[j][k]}")
                    parFile.write('\n')
                parFile.close()

            #update step and tNow and pLuxD
            step+=1

    def analysis_plots(self, subVal = 4, yticInt = 4):
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
            fPath = self.ahl12Path/fName
            data_df = pd.read_csv(str(fPath),header=0, index_col = 0)

            vals = data_df.to_numpy(dtype='float')
            vals_sub = vals[::subVal,:]

            # get an array of mean values at the center (10% of the total region)
            meanCent = np.mean(vals[:,self.mid0:self.mid1], axis = 1)

            # add timeseries to the dict
            timeSerDict[vName]= meanCent

            # add final mCherry conc. to a dict to make the dose response curve
            if vName =="mCherry" and self.tAx[-1]>=360:
                self.doseRes[ahl12] =meanCent[int(360/self.dt)]

            if vName =="mCherry" and self.tAx[-1]<360:
                self.doseRes[ahl12] =meanCent[-1]

        # save time series of all variable in the respective ahl dir
        timeSerDF = pd.DataFrame.from_dict(timeSerDict)
        timeSerDF.to_csv(self.ahl12Path/'timeSerAll.csv')

    def timeSerPlots(self, ):
        """
        to make time-series plots of the variable output
        """
        # make a color dict to assign colors to typeId
        viridis = matplotlib.colormaps.get_cmap('magma')
        cmap0 = viridis(np.linspace(0, 1, 1+len(self.ahl12List)))

        for vName in self.vNames:
            # load mCherry timeseries at all ahl treatments
            normSeries = []
            for ahl12 in self.ahl12List:
                fName = self.resultPath/f"ahl12_{ahl12:.3f}"/'timeSerAll.csv'
                data_df = pd.read_csv(str(fName),header=0, index_col = 0)
                vals = np.array(data_df[vName])
                normSeries.append(vals)

            normArr = np.array(normSeries)

            # make plots without scaling the data to maximum of 1
            fig, ax  = plt.subplots(nrows = 1, ncols = 1, figsize= (2.4,1.5))

            for i, ahl12 in enumerate(self.ahl12List):
                ax.plot(self.tAx, normArr[i,:], linewidth = 1, c=cmap0[i], label = f'{ahl12:.0e} nM')

            ax.set_xlabel("Time (min)")
            ax.set_ylabel("Output (A.U.)")
            if np.amax(normArr)>0:
                ax.set_ylim(0.0,1.02*np.amax(normArr))
            ax.set_xlim(0.0,1.02*self.tAx[-1])
            ax.tick_params(axis='both', which='major', labelsize=6)
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title = r"AHL$_{12}$")

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
            if np.amax(normArr)>0:
                normArr0 = normArr/np.amax(normArr)
            else:
                normArr0 = normArr

            # to make time-series plots with mean values from domain center
            fig, ax  = plt.subplots(nrows = 1, ncols = 1, figsize= (2.4,1.5))

            for i, ahl12 in enumerate(self.ahl12List):
                ax.plot(self.tAx, normArr0[i,:], linewidth = 1, c=cmap0[i], label = f'{ahl12:.0e} nM')

            ax.set_xlabel("Time (min)")
            ax.set_ylabel("Output (A.U.)")
            ax.set_ylim(0.0,1.02)
            ax.set_xlim(0.0,1.02*self.tAx[-1])
            ax.tick_params(axis='both', which='major', labelsize=6)
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title = r"AHL$_{12}$")

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

    def saveDoseResponse(self,):
        # make data array
        doseRes = np.array(list(self.doseRes.items()))

        #save doseRes dict as a pandas DataFrame
        doseDf = pd.DataFrame(self.doseRes.items(), columns = ['ahl12', 'pLas'])
        doseDf.to_csv(self.resultPath/"doseRes.csv")

    def plotDoseResponse(self,AHL6List):
        # to make dose-response plots with mean values from domain center
        fig, ax  = plt.subplots(nrows = 1, ncols = 1, figsize= (1.5,1.5))

        # read data from each AHL6 subdir
        dataDict = {}
        colDict = {}
        colDict[AHL6List[0]] = 'k'
        colDict[AHL6List[1]] = 'r'
        yRes = []
        for AHL6 in AHL6List:
            dirName = f'AHL6_{int(AHL6)}'
            dataDict[AHL6] = pd.read_csv(self.resultPath/dirName/"doseRes.csv", header=0, index_col = 0)
            doseRes = np.array(dataDict[AHL6])
            xAra = doseRes[1:,0]
            yRes.append(doseRes[1:,1])

        # normalize pLas output to 1 together for both AHL6 treatements
        yRes = np.array(yRes)
        yResNorm = yRes/np.amax(yRes)

        for i,AHL6 in enumerate(AHL6List):
            ax.scatter(xAra, yResNorm[i], marker= "o", c=colDict[AHL6], s = 2)
            # ax.plot(xAra, yResNorm[i], colDict[AHL6]+'--', linewidth = 1, label = f'{AHL6:.1f}')

            ####
            # dict to store fit parameters
            fitDict = {"K":[], "n":[], "err_K":[], "err_n":[], "ymin":[], "ymax":[]}

            #fit data to sigmoidFunc - 1) use the x and yAvg to calculate the fit parameters; 2) define to new xFit array and calculate yFit using the
            # obtained fit parameters- opt_K and opt_n
            bnds= ([10**-1, 0.5],[10**3, 3.0])
            mthd = 'trf'

            xdata=np.zeros((len(xAra), 3))
            xdata[:,0]= xAra
            xdata[:,1]= np.amin(yResNorm[i])
            xdata[:,2]= np.amax(yResNorm[i])

            popt, pcov = curve_fit(self.sigmoidFunc, xdata, yResNorm[i], bounds=bnds,  method=mthd)
            perr = np.sqrt(np.diag(pcov))

            xFit = np.logspace(-2,3,50)
            opt_K, opt_n = popt

            fitDict["K"].append(opt_K)
            fitDict["n"].append(opt_n)
            fitDict["err_K"].append(perr[0])
            fitDict["err_n"].append(perr[1])
            fitDict["ymin"].append(np.amin(yResNorm[i]))
            fitDict["ymax"].append(np.amax(yResNorm[i]))

            xdata=np.zeros((len(xFit), 3))
            xdata[:,0]= xFit
            xdata[:,1]= np.amin(yResNorm[i])
            xdata[:,2]= np.amax(yResNorm[i])

            yFit = self.sigmoidFunc(xdata,opt_K,opt_n)

            ax.plot(xFit,yFit, linewidth = 1, c = colDict[AHL6], label = f' {AHL6:.1f}')

            fit_df = pd.DataFrame(fitDict,index=[360])
            fName = f"fitParams_{int(AHL6)}.csv"
            fit_df.to_csv(self.resultPath/fName)

        ax.set_xscale("log")
        ax.set_xlabel(r"AHL$_{12}$ (nM)")
        ax.set_ylabel("Output (A.U.)")
        ax.set_ylim(0.0,1.02)
        ax.set_xlim(xAra[0],1.2*xAra[-1])
        ax.legend(title = r"AHL$_6$ (nM)")

        #remove top and right spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(False)
        plt.tight_layout()

        figName_0 = "ahl12_doseRe.png"
        figName_1 = "ahl12_doseRe.pdf"

        # plt.show()
        plt.savefig(self.resultPath/figName_0, transparent=True, dpi = 300)
        plt.savefig(self.resultPath/figName_1, transparent=True)
        plt.close()

    def sigmoidFunc(self, xdata, K, n):
        x = xdata[:,0]
        ymin=xdata[:,1]
        ymax = xdata[:,2]
        y = ymin + (ymax-ymin)*x**n/(x**n+K**n)
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
    datetime_str=(now.strftime("%Y%m%d_"))

    #make output directory
    dirName=str(datetime_str+ "mod4_NegLoopLasR") # the positive loop model

    resultPath=cwdPath/'data_1D'/dirName
    resultPath.mkdir(mode=0o777, parents=True, exist_ok=True)
    print(f"Created result directory {dirName} at {time.time() - t1} sec ...")

    #define global parameters
    steps = 720
    dt= 1.0 # time in min
    gridSize = 2 # mm
    dx= 0.2 # mm
    nx = int(gridSize/dx)
    vNames = ['LasR', 'mCherry', 'Aqs1']

    # define diffusion rates (per min)
    D6= 0. #3.6*1E-2 # AHL6 diffusion rate(act)
    D12 = 0. #1.8*1E-2 # AHL12 diffusion rate (inh)

    # synthesis rates (1/min, )
    k_prot = 0.2 # approx.
    k_mChr = 0.02

    # degradation rates (unit = 1/min, )
    # note: half-life of 70 min --> degradation rate = 0.01 1/min
    kd_ssrA = 0.018 # half life of 38.5 min - based on http://parts.igem.org/Part:BBa_K1399004
    kd_prot = 0.0028 # half-lie ~ 4 h, Aqs1 doesn't have a ssrA tag

    # activation threshold
    K_T6 = 10. # nM
    K_T12 = 50.0 # nM for LasR-AHL12
    K_T12xR = 200.0 # nM - for LuxR-AHL12
    K_Taqs = 200.0 # nM - for Aqs1-LasR
    nH_12 = 2.0 # Hill coefficient for AHL12
    nH_6 = 1.5 # Hill-C for AHL6
    nH_aqs = 2.0 # Hill C for Aqs1-LasR binding

    k_leak = 0.08 # leaky expression in general

    LuxR = 10.0 # nM
    AHL6List = [0.0, 1000.0] # nM
    ahl12List = [0.0, 0.1, 1., 10., 100., 1E3] #[0.1, 1.0, 10.0, 100., 1000.0] #

    # ahl12 = 1.0
    for AHL6 in AHL6List:
        mod_negL = RD_simulate(nx, dx, dt, steps, vNames, resultPath, ahl12List)
        dirName = f'AHL6_{int(AHL6)}'
        mod_negL.resultPath = resultPath/dirName # to make a subdir for the AHL6 treatment
        mod_negL.resultPath.mkdir(mode=0o777, parents=True, exist_ok=True)

        y0 = [0.1,0.1,0.1] # initial conditions  ['LasR', 'mCherry', 'Aqs1']
        for ahl12 in ahl12List:
            rates = [k_prot, kd_prot, kd_ssrA, k_mChr, K_T6, K_T12, \
            K_T12xR, K_Taqs, nH_6, nH_12, nH_aqs, LuxR, ahl12, AHL6, k_leak]

            mod_negL.makeModel(y0,rates,D6,D12)
            mod_negL.simulateModel()
            mod_negL.analysis_plots(subVal = 4, yticInt = 4)

            print(f'Completed ahl12 dose = {ahl12} at AHL6 = {AHL6}')
        mod_negL.saveDoseResponse()
        mod_negL.timeSerPlots()

    # make a single plot of AHL12 dose-response under AHL6-treated and untreated conditions
    mod_negL.resultPath = resultPath
    mod_negL.plotDoseResponse(AHL6List)

    #end code
    t2= time.time()
    totalSec= t2-t1
    Sec=int(totalSec%60)
    Hrs=int(totalSec//3600)
    Min=int((totalSec%3600)//60)

    print ("Program completed in %sHr:%sMin:%ssec\n"%(Hrs,Min,Sec))
