"""
Amit Landge
2023-07-25

Goal -
To numerically simulate the behaviour of the AHL12 sensor circuit, implmented in E. coli.
In E. coli this circuit has produce AHL12-dose dependent output of fluorescence reporter GFP.

"""

# import necessary packages
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import time
import argparse
import os
from pathlib import Path
from datetime import datetime
import copy
from fipy import *
import pandas as pd
from matplotlib import cm
from scipy.optimize import curve_fit
import seaborn as sns


# define functions
class RD_simulate(object):
    """
    This class contains methods to
    1. initiate the parameters for the reaction-diffusion system.
    2. define the reaction-diffusion PDEs
    3. numerically solve the PDEs
    4. analyze the results - make time-series plots, space-time plots, dose-response plots.
    """

    def __init__(self, nx, dx, dt, steps, vNames, resultPath, ahl12List):
        """
        The init method.
        """
        super(RD_simulate, self).__init__()
        self.mesh= Grid1D(nx = nx, dx = dx)
        self.nx = nx
        self.dx = dx
        self.steps = steps
        self.dt = dt
        self.tAx = [i*dt for i in range(steps+2)]
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

    def makeModel(self,y0,rates,D12):
        """
        method to define PDEs, initial and boundary conditions
        """
        #unpack rates
        k_mChr, kd_mChr, k_lasR, kd_lasR, K_T, nH, k_leak, ahl12  = rates

        # make a separate directory for each arabinose conc results
        self.ahl12Path=self.resultPath/f"ahl12_{ahl12:.3f}"
        self.ahl12Path.mkdir(mode=0o777, parents=True, exist_ok=True)

        # save params to a text file
        with open(self.ahl12Path/"params.csv", "w") as parFile:
            #write the header
            parFile.write("k_mChr, kd_mChr, k_lasR, kd_lasR, K_T, nH, k_leak, ahl12, D12")
            parFile.write('\n')
            parFile.write(f"{k_mChr}, {kd_mChr}, {k_lasR}, {kd_lasR}, {K_T}, {nH}, {k_leak}, {ahl12}, {D12}")
            parFile.write('\n')
        parFile.close()

        # starting state
        self.InVal= y0

        # create CellVariables
        self.vL = [] # list of CellVariables
        for i, vName in enumerate(self.vNames):
            # to add random initial noise
            var = CellVariable(name = vName, mesh = self.mesh, value = self.InVal[i]) # domeShp*
            self.vL.append(var)

        # set initial AHL12
        self.vL[2].setValue(ahl12)

        #define the reaction part
        # note [0 - mCherry, 1 - LasR, 2 - AHL12]
        self.LasR_AHL12 = self.hillF(self.vL[1], self.vL[2], K_T, nH, k_leak)

        eq00 = k_mChr*self.LasR_AHL12- kd_mChr*self.vL[0] # mCherry
        eq01 = k_lasR  - kd_lasR*self.vL[1] # LasR
        eq02 = - kd_ahl*self.vL[2] # AHL12

        # define reaction-diffusion system of PDEs
        eq_mCherry = TransientTerm(var=self.vL[0]) == DiffusionTerm(coeff=0.0, var=self.vL[0])+ eq00
        eq_LasR = TransientTerm(var=self.vL[1]) == DiffusionTerm(coeff=0.0, var=self.vL[1])+ eq01
        eq_AHL12 = TransientTerm(var=self.vL[2]) == DiffusionTerm(coeff=D12, var=self.vL[2])+ eq02

        self.eq = eq_mCherry & eq_LasR & eq_AHL12

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
        method to numerically solve the system and save result to a dictionary
        """
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
            fPath = self.ahl12Path/fName
            data_df = pd.read_csv(str(fPath),header=0, index_col = 0)

            vals = data_df.to_numpy(dtype='float')
            vals_sub = vals[::subVal,:]

            # get an array of mean values at the center (10% of the total region)
            meanCent = np.mean(vals[:,int(0.45*self.nx):int(0.55*self.nx)], axis = 1)

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
        viridis = matplotlib.colormaps.get_cmap('magma', 256)
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
            normArr0 = normArr/np.amax(normArr)

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

    def makeDoseResponse(self,):
        """
        method to make dose-response plots with mean values from domain center
        """
        fig, ax  = plt.subplots(nrows = 1, ncols = 1, figsize= (1.5,1.5))

        # make data array
        doseRes = np.array(list(self.doseRes.items()))

        #save doseRes dict as a pandas DataFrame
        doseDf = pd.DataFrame(self.doseRes.items(), columns = ['AHL12', 'mCherry'])
        doseDf.to_csv(self.resultPath/"doseRes.csv")

        xAra = doseRes[1:,0]
        yRes = doseRes[1:,1]/np.amax(doseRes[:,1])
        ax.scatter(xAra, yRes, marker= "o", c='k', s = 5)

        ####
        # dict to store fit parameters
        fitDict = {"K":[], "n":[], "err_K":[], "err_n":[], "ymin":[], "ymax":[]}

        #fit data to sigmoidFunc - 1) use the x and yAvg to calculate the fit parameters; 2) define to new xFit array and calculate yFit using the
        # obtained fit parameters- opt_K and opt_n
        bnds= ([10**-2, 0.5],[10**2, 2.0])
        mthd = 'trf'

        xdata=np.zeros((len(xAra), 3))
        xdata[:,0]= xAra
        xdata[:,1]= np.amin(yRes)
        xdata[:,2]= np.amax(yRes)

        popt, pcov = curve_fit(self.sigmoidFunc, xdata, yRes, bounds=bnds,  method=mthd)
        perr = np.sqrt(np.diag(pcov))

        xFit = np.logspace(-3,3,50)
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

        yFit = self.sigmoidFunc(xdata,opt_K,opt_n)

        ax.plot(xFit,yFit, c="k", linewidth = 1)
        ####
        # ax.plot(xAra, yRes, 'g--', linewidth = 1)

        ax.set_xscale("log")
        ax.set_xlabel(r"AHL$_{12}$ (nM)")
        ax.set_ylabel("Output (A.U.)")
        ax.set_ylim(0.0,1.02*np.amax(yRes))
        ax.set_xlim(xAra[0],1.2*xAra[-1])

        #remove top and right spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(False)
        plt.tight_layout()

        figName_0 = f"ahl12_doseRe.png"
        figName_1 = f"ahl12_doseRe.pdf"

        # plt.show()
        plt.savefig(self.resultPath/figName_0, transparent=True, dpi = 300)
        plt.savefig(self.resultPath/figName_1, transparent=True)
        plt.close()

        fit_df = pd.DataFrame(fitDict,index=[360])
        fit_df.to_csv(self.resultPath/f"fitParams.csv")

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
    datetime_str=(now.strftime("%Y%m%d_")) # %H%M%S_

    #make output directory
    dirName=str(datetime_str+ "mod2_AHL12sensor_doseRe") # the positive loop model

    resultPath=cwdPath/'data_1D'/dirName
    resultPath.mkdir(mode=0o777, parents=True, exist_ok=True)
    print(f"Created result directory {dirName} at {time.time() - t1} sec ...")

    #define global parameters
    steps = 180
    dt= 2.0 # time in min
    gridSize = 1 # mm
    dx= 0.1 # mm
    nx = int(gridSize/dx)
    vNames = ['mCherry', 'LasR', 'AHL_12']

    """
    *Note: to implement delay in observation of mCherry signal (maturation time ~ 60 min) -
    https://pubs.acs.org/doi/10.1021/acssynbio.1c00387
    I'll simply make another variable mCherry that is closer to the experimental mCherry measurements
    """

    # define diffusion rates (per min)
    D12= 0.0 #3.6*1E-2 # AHL12 diffusion rate(act)

    # synthesis rates (1/min, )
    k_mChr = 0.02 # synthesis time  ~ 5 min
    k_lasR = 0.2 # synthesis time ~ 5 min

    # degradation rates (unit = 1/min, )
    # note: half-life of 60 min --> degradation rate = 0.01 1/min
    kd_mChr = 0.018 # half life of 38.5 min - based on http://parts.igem.org/Part:BBa_K1399004
    # original paper - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC106306/

    kd_lasR = 0.0028 #t1/2 ~ 240 min
    kd_ahl = 0.0011 #t1/2 ~ 10 h = 600 min - based on
    # table 1 of --> https://academic.oup.com/femsec/article/52/1/13/482593

    # activation threshold
    K_T = 50.0 # nM
    nH = 2.0 # Hill coefficient
    k_leak = 0.05 # leaky expression of mCherry from mCherry promoter

    ahl12List = [0.0, 0.01, 0.1, 1., 10., 100., 1000.]
    # ahl12List = [10., 100., 1000.]

    mod_L12 = RD_simulate(nx, dx, dt, steps, vNames, resultPath, ahl12List)

    # y0 = [1.0 for i in range(2)]
    y0 = [1.0, 5., 0.01]

    for ahl12 in ahl12List:
        rates = [k_mChr, kd_mChr, k_lasR, kd_lasR, K_T, nH, k_leak, ahl12]
        mod_L12.makeModel(y0,rates,D12)
        mod_L12.simulateModel()
        mod_L12.analysis_plots()
        print(f'Completed ahl12 dose = {ahl12}')

    mod_L12.makeDoseResponse()
    mod_L12.timeSerPlots()

    #end code
    t2= time.time()
    totalSec= t2-t1
    Sec=int(totalSec%60)
    Hrs=int(totalSec//3600)
    Min=int((totalSec%3600)//60)

    print ("Program completed in %sHr:%sMin:%ssec\n"%(Hrs,Min,Sec))
