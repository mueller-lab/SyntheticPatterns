"""
Amit Landge
2023-07-23

Goal -
To numerically simulate the behaviour of the AHL6 sensor circuit, implmented in E. coli.
In E. coli this circuit has produce AHL6-dose dependent output of fluorescence reporter GFP.

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
# from scipy.optimize import curve_fit
from scipy.optimize import minimize
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

    def __init__(self, nx, dx, dt, steps, vNames, resultPath):
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

    def hillF(self,x,y,Kt, n = 2, k_leak = 0.01):
        """
        Hill function
        """
        return x*(k_leak + y**n/(Kt**n + y**n))

    def getLawnRD(self, secR = 0.5, shift= 0.0):
        """
        a method to define lawn-occupied area of the agar
        reaction-part of the system happens only on the lawn occupied area
        diffusion happens throught entire agar area

        secR: the fraction of the domain length that is occupied by cells
        shift: to shift the lawn-center to right (+) or left (-) as % of domain-length
        """
        secD = 1. - secR # agar-area fraction without the cell-lawn

        leftInd = int((0.5*secD+shift)*self.nx)
        rightInd = int(((1-0.5*secD)+shift)*self.nx)
        LawnArea = np.zeros(self.nx)
        LawnArea[leftInd:rightInd] = 1.0
        # LawnArea[:] = 1.0
        self.LawnArea = LawnArea.astype('float')

        # define lawn center for inducer addition
        centL = int(0.44*self.nx)
        centR = int(0.56*self.nx)
        self.LawnCent= np.zeros(self.nx)
        self.LawnCent[centL:centR] = 1.

        # make LawnCent a domeShape curve instead of the step-function
        centLen = int(np.sum(self.LawnCent))
        # angleRang = np.linspace(0., np.pi,centLen)
        # domeShape= np.sin(angleRang)
        linShape1 = np.linspace(1., 0., int(centLen/2))
        linShape2 = np.linspace(0., 1., centLen - int(centLen/2))

        linShape = np.concatenate((linShape2, linShape1))
        self.LawnCent[centL:centR] = linShape
        # self.LawnCent[centL:centR] = domeShape

        # to make initial dome-shaped trigger
        domeShp = np.zeros(self.nx)
        domeShape = np.linspace(0.0, np.pi, int(secR*self.nx))
        domeShape= np.sin(domeShape) # make a domeShape curve using sine wave

        # add the dome-shpae to the initial domeShp array - make sure domeShape has same length as the section of domeShp
        lendS = len(domeShape)
        domeShp[leftInd: lendS + leftInd] = domeShape
        self.domeShp = domeShp.astype('float')

    def makeModel(self,y0,rates,secR, D6):
        """
        method to define PDEs, initial and boundary conditions
        """
        #unpack rates
        k_gfp, kd_ssrA, k_luxR, kd_luxR, K_T, nH, k_leak, ahl6  = rates

        # make a separate directory for each arabinose conc results
        self.ahl6Path=self.resultPath/f"ahl6_{ahl6:.3f}"
        self.ahl6Path.mkdir(mode=0o777, parents=True, exist_ok=True)

        # save params to a text file
        with open(self.ahl6Path/"params.csv", "w") as parFile:
            #write the header
            parFile.write("k_gfp, kd_ssrA, k_luxR, kd_luxR, K_T, nH, k_leak, ahl6, D6")
            parFile.write('\n')
            parFile.write(f"{k_gfp}, {kd_ssrA}, {k_luxR}, {kd_luxR}, {K_T}, {nH}, {k_leak}, {ahl6}, {D6}")
            parFile.write('\n')
        parFile.close()

        # starting state
        self.InVal= y0

        #get LawnArea
        self.getLawnRD(secR)

        # create CellVariables
        self.vL = [] # list of CellVariables
        for i, vName in enumerate(self.vNames):
            # to add random initial noise
            var = CellVariable(name = vName, mesh = self.mesh, value = self.InVal[i]*self.LawnArea) # domeShp*
            self.vL.append(var)

        # set initial AHL6
        self.vL[2].setValue(ahl6*self.LawnCent)

        #define the reaction part
        # note [0 - sfGFP, 1 - LuxR, 2 - AHL6]
        self.LuxR_AHL6 = self.hillF(self.vL[1], self.vL[2], K_T, nH, k_leak)

        eq00 = self.LawnArea*(k_gfp*self.LuxR_AHL6- kd_ssrA*self.vL[0]) # sfGFP
        eq01 = self.LawnArea*(k_luxR  - kd_luxR*self.vL[1]) # LuxR
        eq02 = - kd_ahl*self.vL[2] # AHL6

        # define reaction-diffusion system of PDEs
        eq_sfGFP = TransientTerm(var=self.vL[0]) == DiffusionTerm(coeff=0.0, var=self.vL[0])+ eq00
        eq_LuxR = TransientTerm(var=self.vL[1]) == DiffusionTerm(coeff=0.0, var=self.vL[1])+ eq01
        eq_AHL6 = TransientTerm(var=self.vL[2]) == DiffusionTerm(coeff=D6, var=self.vL[2])+ eq02

        self.eq = eq_sfGFP & eq_LuxR & eq_AHL6

        # #write results to .csv files
        # for i, vName in enumerate(self.vNames):
        #     fName = f"{vName}.csv"
        #     with open(self.ahl6Path/fName, "w") as parFile:
        #         #write the header
        #         parFile.write("t_min")
        #         for j in range(self.nx):
        #             parFile.write(f",{j*self.dx}")
        #         parFile.write('\n')
        #
        #         #write the first row
        #         parFile.write("0.00")
        #         for j in range(self.nx):
        #             parFile.write(f",{self.vL[i][j]}")
        #         parFile.write('\n')
        #     parFile.close()

    def simulateModel(self, ):
        """
        method to numerically solve the system and save result to a dictionary
        """
        step = 0
        self.fluorDict = {}

        while step <= self.steps:
            # add fluorescent readOut at selected timepoints to a dictionary
            # 0, 1, 2, 3, 4 h
            tNow = step*self.dt
            if (tNow%60)==0:
                hr = int(tNow/60)
                self.fluorDict[hr] = np.array(self.vL[0])

            self.eq.solve(dt=self.dt)

            # #save results to the .csv files
            # for j, vName in enumerate(self.vNames):
            #     fName = f"{vName}.csv"
            #     with open(self.ahl6Path/fName, "a") as parFile:
            #         #append a new row to the file
            #         parFile.write(f"{self.dt*(step+1):.2f}")
            #         for k in range(self.nx):
            #             parFile.write(f",{self.vL[j][k]}")
            #         parFile.write('\n')
            #     parFile.close()
            step+=1

        return self.fluorDict

    def plot_fluor(self,):
        """
        method to -
        1. plot fluorescent output every hour
        2. save fluorDict to csv
        """
        fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize= (1.6,1.6))
        x = np.linspace(0, self.nx*self.dx, self.nx)# length in mm

        # make a color dict to assign colors to timepoints
        viridis = cm.get_cmap('viridis', 256)
        hrCol = viridis(np.linspace(0, 1, len(self.fluorDict.keys())))

        #get global max of all values in the dict
        globMax = np.max(np.array(list(self.fluorDict.values())))

        for key, val in self.fluorDict.items():
            ax.plot(x, val/globMax, linewidth=1, c=hrCol[key], label = f"{key} h")

        ax.set_xlabel("r (mm)")
        ax.set_ylabel("Output (A.U.)")
        ax.set_ylim(0,1.02)
        ax.set_xlim(0, x[-1]*1.02)
        ax.tick_params(axis='both', which='major', labelsize=6)
        ax.legend()
        plt.legend(frameon=False)

        #remove top and right spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        plt.tight_layout()
        figName_0 = "fluor_profile.png"
        figName_1 = "fluor_profile.pdf"
        # plt.show()
        plt.savefig(self.resultPath/figName_0, transparent=True, dpi = 300.)
        plt.savefig(self.resultPath/figName_1, transparent=True)
        plt.close()

        # save data to csv
        fluorDf = pd.DataFrame(self.fluorDict)
        fluorDf.to_csv(self.resultPath/'fluor_profile.csv')

    def plotRadProf(self, xAx, dataArr, headIndex = 0, makePlot= False, hrRange = range(5)):
        """
        method to -
        1. return radial profile - comparable to the experimetal data
        2. if makePlot==True, plot fluorescent output radial profile every hour
        """
        xCent = int(0.5*self.nx)+headIndex # center index for the lawn
        # x = np.linspace(0, xLen*self.dx, xLen)# length in mm
        xLen = len(xAx)

        # #get global max of all values in the dict
        # globMax = np.max(np.array(list(self.fluorDict.values())))

        simDataArr = []
        for key, val in self.fluorDict.items():
            simDataArr.append(val[xCent:xLen+xCent])

        simDataArr= np.array(simDataArr)

        # select data in the hrRange
        hrList = [i for i in hrRange]
        simDataArr = simDataArr[hrList]
        dataArr = dataArr[hrList]

        # normalize at each timePoint
        for i in range(dataArr.shape[0]):
            simDataArr[i] = simDataArr[i]/np.max(simDataArr[i])
            dataArr[i] = dataArr[i]/np.max(dataArr[i])

        # make a color dict to assign colors to timepoints
        if makePlot:
            fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize= (2.4,1.6))
            viridis = cm.get_cmap('viridis', 256)
            hrCol = viridis(np.linspace(0, 1, len(hrList)))

            for i, val in enumerate(simDataArr):
                ax.plot(xAx, dataArr[i], linewidth=1, c=hrCol[i], label = f"{hrList[i]} h") # plot experimetal
                ax.plot(xAx, val, '--', linewidth=1, c=hrCol[i], label = f"{hrList[i]} h [fit]", alpha = 0.5) # plot simulated

            ax.set_xlabel("r (mm)")
            ax.set_ylabel("Output (A.U.)")
            ax.set_ylim(0,1.02)
            ax.set_xlim(0, xAx[-1]*1.02)
            ax.tick_params(axis='both', which='major', labelsize=6)
            plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

            #remove top and right spines
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

            plt.tight_layout()
            figName_0 = "sim_fluor_profile.png"
            figName_1 = "sim_fluor_profile.pdf"
            # plt.show()
            plt.savefig(self.resultPath/figName_0, transparent=True, dpi = 300.)
            plt.savefig(self.resultPath/figName_1, transparent=True)
            plt.close()

        return simDataArr, dataArr

def getExptProf(drop, numHr, outDir, fluor = 'gfp', chopOff=200.):
    """
    method to read the experimetal radial profile data and
    plot the normalized data at all the time points
    """
    dataDict={}
    drop= str(drop).zfill(2)

    for i in range(numHr):
        fName = "%s_%sh-%s.csv"%(fluor, i, drop)
        pdDf = pd.read_csv(dataDir/fName, header=0, index_col = 0)
        dataDict[i] = pdDf

    # processing data
    xAx = np.array(dataDict[0].index)
    # print(xAx.shape)

    # remove the head of the gradient where the drop is added -
    # this head part gives a noisy fluorescence readout
    # From the known drop size it looks that 200 um part should be removed
    headIndex= np.where(chopOff<xAx)[0][0]

    xAx = xAx[headIndex:]-xAx[headIndex] # now the xAx is shifted to left by 200 um
    # print(xAx.shape)
    # print(xAx[0])

    # now remove first (headIndex, eg. 149) values from the dataDict for all the time points
    for x, y in dataDict.items():
        dataDict[x] = y.iloc[headIndex:,0]

    # subtract autofluorescence (the minima of the 0 h readings)
    real_0h = dataDict[0]
    min_0h = np.min(real_0h)
    for x, y in dataDict.items():
        dataDict[x]= y-min_0h

    # normalize the data by dividing by max value
    globMax = np.max(np.array(list(dataDict.values()))) # get global max
    normDict = {}
    for x, y in dataDict.items():
        normDict[x] = y/globMax # global normalization

    # convert data to numpy.array format
    dataArr = np.array(list(normDict.values()))

    # subtract autofluorescence (0 h data) from dataArr
    for i in range(dataArr.shape[0]):
        dataArr[i] = dataArr[i] - dataArr[0]

    # plot normalized data at different time-points
    # make figure
    fig, ax  = plt.subplots(nrows = 1, ncols = 1, figsize= (1.6,1.6))

    # make a color dict to assign colors to typeId
    viridis = cm.get_cmap('viridis', 256)
    cmap0 = viridis(np.linspace(0, 1, numHr))

    xAx_mm = xAx/1000 # xAx values in mm

    # j = 0
    for i in range(0,numHr):
        ax.plot(xAx_mm, dataArr[i], label= "%s h"%(i), c = cmap0[i], linewidth = 1.)

    plt.xlabel("Length (mm)")
    plt.ylabel("Normalized intensity (A.U.)")
    plt.legend(frameon=False)

    #set axis limits
    ax.set_ylim(0.,1.02)
    ax.set_xlim(0., xAx_mm[-1]*1.02)

    ax.grid(False)
    plt.tight_layout()

    #remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    fName = "%s_drop%s.pdf"%(fluor, drop)
    fName1 = "%s_drop%s.png"%(fluor, drop)

    plt.savefig(outDir/fName, transparent= True)
    plt.savefig(outDir/fName1, transparent= True, dpi = 300.)
    plt.close()
    # plt.show()
    return dataArr, xAx_mm, headIndex

def objFun(x, dataArr=None, xAx=None, headIndex=210, hrRange= range(5), makePlot= False):
    """
    cost function or objective function to use for minimization of the residual using
    scipy.optimize.minimize
    """
    #run simulation
    mod_L6 = RD_simulate(nx, dx, dt, steps, vNames, resultPath) #
    mod_L6.makeModel(y0,rates,secR,x) # x = D to be optimized
    fluorDict = mod_L6.simulateModel()
    simDataArr, dataArr = mod_L6.plotRadProf(xAx, dataArr,headIndex,makePlot, hrRange)

    # calculate RMSD
    MSDiff = [(dataArr[i] - simDataArr[i])**2 for i in range(len(simDataArr))]
    tot_MSDiff=[np.sum(j) for j in MSDiff]
    tot_Diff = np.sqrt(sum(tot_MSDiff)/len(tot_MSDiff))
    print(f'tot_Diff = {tot_Diff:.2f}')
    return tot_Diff

# main code
if __name__=="__main__":
    # get current dir and create result dir
    print("\n**************************************************** \nStarting the program...")
    t1 = time.time()

    # global plot settings
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["font.size"] = 6
    plt.rcParams['legend.fontsize'] = 6

    #get current dir and datetime
    cwdPath=Path(os.path.abspath(os.getcwd()))
    now=datetime.now()
    datetime_str=(now.strftime("%Y%m%d_")) # %H%M%S_

    #define global parameters
    steps = 60
    dt= 3.0 # time in min
    gridSize = 30 # mm
    dx= 19.017/2048 # mm 19.017/2048 gives the experimetal pixel size in mm
    nx = int(gridSize/dx)
    vNames = ['sfGFP', 'LuxR', 'AHL_6']


    # define diffusion rates (per min)
    D6= 3.0*1E-2 # AHL6 diffusion rate(act)

    # synthesis rates (1/min, )
    k_gfp = 0.1 # synthesis time  ~ 5 min
    k_luxR = 0.2 # synthesis time ~ 5 min

    # degradation rates (unit = 1/min, )
    # note: half-life of 60 min --> degradation rate = 0.01 1/min
    kd_ssrA = 0.018 # half life of 38.5 min - based on http://parts.igem.org/Part:BBa_K1399004
    # original paper - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC106306/

    kd_luxR = 0.0028 #t1/2 ~ 240 min
    kd_ahl = 0.0011 #t1/2 ~ 10 h = 600 min - based on
    # kd_ahl = kd_ahl+0.002 # to account for out-of-plane loss of ahls in the agar
    # table 1 of --> https://academic.oup.com/femsec/article/52/1/13/482593

    # activation threshold
    K_T = 10.0 # nM
    nH = 1.5 # Hill coefficient
    k_leak = 0.05 # leaky expression of sfGFP from sfGFP promoter
    # ahl6 = 400.
    ahl6 = 400. # 2021-03-05

    # y0 = [1.0 for i in range(2)]
    y0 = [0.01,0.01,0.01]
    rates = [k_gfp, kd_ssrA, k_luxR, kd_luxR, K_T, nH, k_leak, ahl6]
    secR = 2/3

    #make output directory
    dropList = [5,6] # 2021-04-16
    # dropList = [1,2,3,4,5,6] # 2021-03-11
    # dropList = [2,3,4,5] # 2021-03-05

    datePath = cwdPath/'out'/'2021-04-16'
    if not(datePath.exists()):
        datePath.mkdir(parents=True)

    with open(datePath/'out.csv', 'w') as f:
        f.write('drop, D (um^2/s), residual \n')

    for drop in dropList:
        dirName=str(datetime_str+ f"mod1_AHL6sensor_drop{drop}_subAuto_Fig3") # the positive loop model

        resultPath=cwdPath/'out'/'2021-04-16'/dirName
        resultPath.mkdir(mode=0o777, parents=True, exist_ok=True)
        print(f"Created result directory {dirName} at {time.time() - t1} sec ...")


        # get the experimetal radial profile
        dataDir = cwdPath/"data"/"2021-04-16" # input csv file path
        # dataDir = cwdPath/"data"/"2021-03-05" # input csv file path
        numHr = 5
        outDir = resultPath
        dataArr, xAx, headIndex = getExptProf(drop, numHr, outDir, fluor = 'gfp', chopOff=2800.)

        print(f'headIndex = {headIndex}')

        # test run the simulation with the experimetal settings
        # mod_L6 = RD_simulate(nx, dx, dt, steps, vNames, resultPath)
        # mod_L6.makeModel(y0,rates,secR,D6)
        # fluorDict = mod_L6.simulateModel()
        # simDataArr = mod_L6.plotRadProf(xAx, dataArr,headIndex, makePlot= True)

        # test run objFun
        hrRange = range(1,4)
        # tot_Diff = objFun(0.022, dataArr, xAx, headIndex, hrRange= hrRange, makePlot = True)

        # use scipy.optimize.minimize to find optimal value of D
        x0 = 500. * 60/1E6
        bnds= ((0.0010,0.050),)
        res = minimize(objFun, x0=x0, args=(dataArr, xAx, headIndex,hrRange), bounds=bnds, \
        options = {'maxiter':20})

        # convert D to um^2/s
        D_um = res.x[0]*1E6/60

        # simulate model with optimized parameters
        tot_Diff = objFun(res.x, dataArr, xAx, headIndex, hrRange= hrRange, makePlot = True)

        with open(datePath/'out.csv', 'a') as f:
            f.write(f'{drop}, {D_um}, {tot_Diff} \n')


    #end code
    t2= time.time()
    totalSec= t2-t1
    Sec=int(totalSec%60)
    Hrs=int(totalSec//3600)
    Min=int((totalSec%3600)//60)

    print ("Program completed in %sHr:%sMin:%ssec\n"%(Hrs,Min,Sec))
