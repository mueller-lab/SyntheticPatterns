"""
Amit Landge

Goal -To model the complete synthetic patterning circuit (on a 2D grid)

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

# define functions
class RD_simulate(object):
    """docstring for RD_simulate."""

    def __init__(self, nx, dx, dt, steps, vNames, resultPath, araList, secR, shift, uniF):
        super(RD_simulate, self).__init__()
        # self.mesh= Grid1D(nx = nx, dx = dx)
        # self.mesh= Grid2D(nx=nx, ny=nx, dx=dx, dy=dx)
        cellSize = dx
        radius = int(nx*dx/2)
        self.mesh = Gmsh2D('''
                        cellSize = %(cellSize)g;
                        radius = %(radius)g;
                        Point(1) = {0, 0, 0, cellSize};
                        Point(2) = {-radius, 0, 0, cellSize};
                        Point(3) = {0, radius, 0, cellSize};
                        Point(4) = {radius, 0, 0, cellSize};
                        Point(5) = {0, -radius, 0, cellSize};
                        Circle(6) = {2, 1, 3};
                        Circle(7) = {3, 1, 4};
                        Circle(8) = {4, 1, 5};
                        Circle(9) = {5, 1, 2};
                        Line Loop(10) = {6, 7, 8, 9};
                        Plane Surface(11) = {10};
                        '''  % locals() )

        self.nx = nx
        self.dx = dx
        self.Xc, self.Yc = self.mesh.cellCenters
        self.Xca = np.array(self.Xc)
        self.Yca = np.array(self.Yc)
        self.cirRds = radius

        self.tAx = [i*dt for i in range(steps+2)]
        self.steps = steps
        self.dt = dt
        self.resultPath = resultPath
        self.vNames = vNames
        self.nodeN = len(vNames)
        self.doseRes = {}
        self.araList = araList
        self.secR = secR
        self.shift = shift
        self.tNow = 0
        self.LasR = 40.
        self.LuxR = 40.
        self.LasI = 2000.0
        self.uniF = uniF

    def hillF(self,x,y,Kt, n = 2, k_leak = 0.01):
        """
        Hill function
        """
        return x*(k_leak + y**n/(Kt**n + y**n))

    def makeCirMask(self,):
        # to restrict the reactions (cells) in the middle circular domain
        maskRad = self.cirRds*self.secR
        maskCent = 0.0 + self.shift
        self.CirMask  = (((self.Xc - maskCent)**2 + (self.Yc - maskCent)**2)<maskRad**2)

        radPosArr= np.sqrt(((self.Xca- maskCent)**2 + (self.Yca- maskCent)**2)/maskRad**2)
        radPosArr[radPosArr > 1.0]= 1.0
        self.iniDome =0.2+ np.sin(0.5*np.pi*(1.0 - radPosArr))
        self.iniDome[~self.CirMask] = 0.0

        if self.uniF:
            self.iniDome[self.iniDome>0] =1.0

    def makeModel(self,y0,ara,D6,D12,Dc):
        # make a separate directory for each ara conc results
        self.araPath=self.resultPath/f"ara_{ara:.2f}"
        self.araPath.mkdir(mode=0o777, parents=True, exist_ok=True)

        # save params to a text file
        with open(self.araPath/"params.csv", "w") as parFile:
            #write the header
            parFile.write("k_prot, k_gfp, k_mChr, kd_prot, kd_ssrA, k_ahl, kd_ahl, kda_ahl6, kda_ahl12, k_cell, K_cell, \
            K_T6, K_T12, K_T12xR, K_Taqs, K_Tara, k_leak, ara, \
            nH_6, nH_12, nH_ara, nH_aqs, secR, shift, LasI, LasR, LuxR")
            parFile.write('\n')
            parFile.write(f"{k_prot}, {k_gfp}, {k_mChr}, {kd_prot}, {kd_ssrA}, {k_ahl}, {kd_ahl}, {kda_ahl6}, {kda_ahl12}, \
            {k_cell}, {K_cell}, {K_T6}, {K_T12}, {K_T12xR}, {K_Taqs}, {K_Tara}, {k_leak}, {ara}, \
            {nH_6}, {nH_12}, {nH_ara}, {nH_aqs}, {self.secR}, {self.shift}, {self.LasI}, {self.LasR}, {self.LuxR}")
            parFile.write('\n')

        # starting state
        self.InVal= y0

        # create CellVariables
        self.vL = [] # list of CellVariables
        for i, vName in enumerate(self.vNames):
            var = CellVariable(name = vName, mesh = self.mesh, value = self.InVal[i]*self.iniDome) #*self.iniDome
            self.vL.append(var)

        #define the reaction part -
        # ['LuxI', 'AHL12', 'AHL6', 'mCherry', 'sfGFP', 'AiiA', 'Aqs1','cellD']
        #   0       1         2       3         4       5         6        7
        self.LasR_inh = self.hillF(1.0, K_Taqs, self.vL[6], nH_aqs, k_leak)
        self.LasR_AHL12 = self.hillF(self.LasR, self.vL[1], K_T12/self.LasR_inh, nH_12, k_leak)
        self.LuxR_AHL6 = self.hillF(self.LuxR,self.vL[2],K_T6,nH_6, k_leak)
        self.LuxR_AHL12 = self.hillF(self.LuxR, self.vL[1], K_T12xR, nH_12, k_leak)

        self.ara_indc = self.hillF(10.,ara,K_Tara,nH_ara,k_leak)

        # self.limitL12 = (1000.0 - self.vL[1])/1000.0
        self.maxLuxI = 800.0
        # self.limitL6 =((self.maxLuxI-self.vL[0])/self.maxLuxI)
        # self.AHL6cap = (K_cell- self.vL[7])/K_cell
        # self.AHL12cap = (1.1*K_cell- self.vL[7])/(1.1*K_cell)
        self.AHL6cap = (K_cell- self.vL[7])/(self.secR*K_cell)
        self.AHL12cap = (1.1*K_cell- self.vL[7])/(1.1*self.secR*K_cell)

        self.pLux = self.LuxR_AHL6 + self.LuxR_AHL12
        self.pLas = self.LasR_AHL12

        self.eq00 = self.vL[7]*(k_prot*(self.pLux + self.pLas) - kd_ssrA*self.vL[0]) # LuxI

        # self.eq01 = self.vL[7]*(k_ahl*self.LasI*self.AHL12cap - self.kd2_ahl*self.vL[5]*self.vL[1]) - kd_ahl*self.vL[1] # AHL12
        # self.eq02 = self.vL[7]*(k_ahl*self.vL[0]*self.AHL6cap - kd1_ahl*self.vL[5]*self.vL[2]) - kd_ahl*self.vL[2] # AHL6

        self.eq01 = self.vL[7]*(k_ahl*self.LasI*self.AHL12cap - kda_ahl12*self.vL[5]*self.vL[1]) - kd_ahl*self.vL[1] # AHL12
        self.eq02 = self.vL[7]*(k_ahl*self.vL[0]*self.AHL6cap - kda_ahl6*self.vL[5]*self.vL[2]) - kd_ahl*self.vL[2] # AHL6

        self.eq03 = self.vL[7]*(k_mChr*self.pLas - kd_ssrA*self.vL[3]) # mCherry
        self.eq04 = self.vL[7]*(k_gfp*self.pLux- kd_ssrA*self.vL[4]) # sfGFP
        self.eq05 = self.vL[7]*(k_prot*self.ara_indc- kd_ssrA*self.vL[5]) # AiiA
        self.eq06 = self.vL[7]*(k_prot*self.pLux- kd_prot*self.vL[6]) # Aqs1
        self.eq07 = self.vL[7]*(k_cell/K_cell)*(K_cell-self.vL[7]) # cellD - logistic cell growth


        # define reaction-diffusion system of PDEs
        self.eq_LuxI = TransientTerm(var=self.vL[0]) == DiffusionTerm(coeff=0.0, var=self.vL[0])+ self.eq00
        self.eq_AHL12 = TransientTerm(var=self.vL[1]) == DiffusionTerm(coeff=D12, var=self.vL[1])+ self.eq01
        self.eq_AHL6 = TransientTerm(var=self.vL[2]) == DiffusionTerm(coeff=D6, var=self.vL[2])+ self.eq02
        self.eq_mChr = TransientTerm(var=self.vL[3]) == DiffusionTerm(coeff=0.0, var=self.vL[3])+ self.eq03
        self.eq_sfGFP = TransientTerm(var=self.vL[4]) == DiffusionTerm(coeff=0.0, var=self.vL[4])+ self.eq04
        self.eq_AiiA = TransientTerm(var=self.vL[5]) == DiffusionTerm(coeff=0.0, var=self.vL[5])+ self.eq05
        self.eq_aqs = TransientTerm(var=self.vL[6]) == DiffusionTerm(coeff=0.0, var=self.vL[6])+ self.eq06
        self.eq_cellD = TransientTerm(var=self.vL[7]) == DiffusionTerm(coeff=Dc, var=self.vL[7])+ self.eq07

        self.eq = self.eq_LuxI & self.eq_AHL12 & self.eq_AHL6 & self.eq_mChr & self.eq_sfGFP & self.eq_AiiA & self.eq_aqs & self.eq_cellD

        #write results to .csv files -
        self.csvPath=self.araPath/"CSVs"
        self.csvPath.mkdir(mode=0o777, parents=True, exist_ok=True)

        # to save min and max values of vNames
        self.minDict = {}
        self.maxDict = {}

        for i, vName in enumerate(self.vNames):
            fName = f"{vName}_0000.csv"
            varArr= np.array(self.vL[i])
            varArr=varArr.astype('float')
            varDf = pd.DataFrame(varArr)
            varDf.to_csv(self.csvPath/fName)

            # save min and max for plotting range
            self.minDict[vName] = 0
            self.maxDict[vName] = np.max(varArr)

    def simulateModel(self, ):
        """
        use eq.solve to iteratively solve the system
        """
        # iterate to numerically solve the system and save result to a dictionary
        step = 0
        self.tNow = 0

        while step <= self.steps:
            self.eq.solve(dt=self.dt)

            # Limit maximum LuxI level to 180
            self.vL[0].setValue(self.maxLuxI, where=self.vL[0]>self.maxLuxI)

            # stop negative conc. in general
            for i in range(8):
                self.vL[i].setValue(0.0, where=self.vL[i]<0.0001)

            #save results to the .csv files
            for i, vName in enumerate(self.vNames):
                stID = f"{step+1}".zfill(4)
                fName = f"{vName}_{stID}.csv"
                varArr= np.array(self.vL[i])
                varArr=varArr.astype('float')

                self.minDict[vName] = min([self.minDict[vName], np.min(varArr)])
                self.maxDict[vName] = max([self.maxDict[vName], np.max(varArr)])

                varDf = pd.DataFrame(varArr)
                varDf.to_csv(self.csvPath/fName)

            #update step and tNow
            step+=1
            self.tNow+=self.dt

            #keep track of progress
            if step%50 ==0:
                print(f'completed step = {step}')

    def analysis_plots(self, ):
        # plot the results as 2D graphs with imShow and save
        tAx = np.arange(0,(self.steps+1)*self.dt,self.dt)
        xAx = np.arange(0,self.nx*self.dx, self.dx)

        # dict to save dose-response data
        self.doseRes[ara] = []

        for j, vName in enumerate(self.vNames):
            #code to plot conc. vs length

            # assign colormap
            if 'mCherry' in vName:
                cmap0 = 'magma'
            else:
                cmap0 = 'viridis'

            #make a new subdir to save plots
            tSeriesPath=self.araPath/f"{vName}_tSeries"
            tSeriesPath.mkdir(mode=0o777, parents=True, exist_ok=True)

            for step in range(0, self.steps+1):
                stID = f"{step}".zfill(4)
                fName = f"{vName}_{stID}.csv"

                fPath = self.csvPath/fName
                data_df = pd.read_csv(str(fPath),header=0, index_col = 0)

                vals = data_df.to_numpy(dtype='float')

                fig, ax  = plt.subplots(nrows = 1, ncols = 1, figsize= (2.2,2.2))
                # pos = ax.imshow(vals, origin = 'lower', cmap =cmap0, vmin = self.minDict[vName], vmax = self.maxDict[vName]) #, vmin = vmin, vmax = vmax
                pos = ax.scatter(self.Xca, self.Yca, c = vals, cmap = cmap0, marker = "h", s = 0.4, edgecolors=None, \
                vmin =self.minDict[vName], vmax =self.maxDict[vName])

                ax.set_xlim([-self.cirRds, self.cirRds])
                ax.set_ylim([-self.cirRds, self.cirRds])
                ax.set_xlabel("Length (mm)", loc = 'center', labelpad = 1.0)
                ax.set_ylabel("Length (mm)", loc = 'center', labelpad = 1.0)
                ax.tick_params(axis='both', which='major', labelsize=6)

                # add colorbar
                # plt.colorbar(pos, fraction = 0.03, pad = 0.2, label = 'Intensity val.', orientation = 'horizontal')
                ax.set_aspect('equal')
                #remove all spines
                ax.spines['top'].set_visible(False)
                # ax.spines['bottom'].set_visible(False)
                ax.spines['right'].set_visible(False)
                # ax.spines['left'].set_visible(False)

                plt.tight_layout()

                figName_0 = f"{vName}_{stID}.png"
                # figName_1 = f"{vName}_{stID}.pdf"

                # plt.show()
                plt.savefig(tSeriesPath/figName_0, transparent=True, dpi = 300)
                # plt.savefig(tSeriesPath/figName_1, transparent=True)
                plt.close()

    def saveDoseResponse(self,):
        # # make data array
        # doseRes = np.array(list(self.doseRes.items()))
        #
        # #save doseRes dict as a pandas DataFrame
        # doseDf = pd.DataFrame(self.doseRes.items(), columns = ['ara', 'pLas'])
        # doseDf.to_csv(self.resultPath/"doseRes.csv")

        # make data array
        xAra = np.array(list(self.doseRes.keys()))
        yArr = np.array(list(self.doseRes.values()))
        # print(yArr)

        #save doseRes dict as a pandas DataFrame
        doseDf = pd.DataFrame(yArr, columns =self.vNames, index = xAra)
        doseDf.to_csv(self.resultPath/"doseRes.csv")

    def plotDoseResponse(self,):
        # to make dose-response plots with mean values from domain center
        fig, ax  = plt.subplots(nrows = 1, ncols = 1, figsize= (2.0,2.0))

        # read data
        dataDF = pd.read_csv(self.resultPath/"doseRes.csv", header=0, index_col = 0)
        doseRes = np.array(dataDF)
        xAra = np.array(dataDF.index)
        yRes = np.array(dataDF['sfGFP_fl'])

        # normalize pLas output to 1
        yResNorm = yRes/np.amax(yRes)

        ax.scatter(xAra, yResNorm, marker= "o", c='k', s = 2)
        ax.plot(xAra, yResNorm, 'g--', linewidth = 1)

        ax.set_xscale("log")
        ax.set_xlabel("Ara (nM)")
        ax.set_ylabel("Output (A.U.)")
        ax.set_ylim(0.0,1.1)
        # ax.set_ylim(0.0,8.0)
        ax.set_xlim(xAra[0],1.1*xAra[-1])

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
    datetime_str=(now.strftime("%Y%m%d_%H%M%S_"))

    #make output directory
    # dirName=str(datetime_str+ "mod5_compCir_circle_Original")
    dirName=str(datetime_str+ "mod5_compCir_circle_Big")

    resultPath=cwdPath/'data_2D'/dirName
    resultPath.mkdir(mode=0o777, parents=True, exist_ok=True)
    print(f"Created result directory {dirName} at {time.time() - t1} sec ...")

    #define global parameters
    steps = 360
    dt= 4.0 # time in min
    gridSize = 16 # mm
    dx= 0.2 # mm
    nx = int(gridSize/dx)
    vNames = ['LuxI', 'AHL12', 'AHL6', 'mCherry', 'sfGFP', 'AiiA', \
    'Aqs1','cellD']

    # define diffusion rates (per min)
    D6= 2.1*1E-2 # AHL6 diffusion rate(act)
    D12 = 2.1*1E-2 # AHL12 diffusion rate (inh)

    # synthesis rates (1/min, )
    k_prot = 0.2 # synthesis time  ~ 5 min
    k_gfp = 0.1 # synthesis + maturation time ~ 10 min
    k_mChr = 0.02 # synthesis + maturation time ~ 50 min
    k_ahl = 0.002 # AHL synthesis rate (by LuxI)

    # degradation rates (unit = 1/min, )
    # note: half-life of 70 min --> degradation rate = 0.01 1/min
    kd_ssrA = 0.018 # half life of 38.5 min - based on http://parts.igem.org/Part:BBa_K1399004
    kd_prot = 0.0056 #t1/2 ~ 2 h
    kd_ahl = 0.002 # ahl decay (hydroslysis, loss into agar), half-life ~ 5 h
    kda_ahl6 = 0.01 # ahl6 degradation by lactonase (arabinose-induced)
    kda_ahl12 = 0.001 # ahl12 degradation by lactonase (arabinose-induced)

    # activation thresholds
    K_T12 = 50.0 # nM for LasR-AHL12
    K_T6 = 10. # nM
    K_T12xR = 2000.0 # nM - for LuxR-AHL12
    K_Taqs = 300.0 # nM - for Aqs1-LasR
    K_Tara = 6E3 # arabinose induction threshold

    nH_12 = 2.0 # Hill coefficient for AHL12
    nH_6 = 1.5 # Hill-C for AHL6
    nH_aqs = 2.0 # Hill C for Aqs1-LasR binding
    nH_ara = 1.5

    k_leak = 0.05 # leaky expression in general

    # to add the effect of cell density
    k_cell = 0.007 # cellD increase per min
    K_cell = 1.6 # cell density maximum capacity
    Dc = 6E-6 # cell diffusion in mm^2/min - calculated from areaIncrease in 26 h- expt. on 2023.01.19

    secR = 1.0 # 0.65 # area-fraction of agar occupied by cell lawn
    shift = 0.0 # a parameter to shift the lawn center ( in mm)
    uniF = False #True
    # uniF = True

    # arabinose treatments
    araList = [0.0] # , 1E6

    mod_compCir = RD_simulate(nx, dx, dt, steps, vNames, resultPath, araList, secR, shift, uniF)
    mod_compCir.makeCirMask()

    # ['LuxI', 'AHL12', 'AHL6', 'mCherry', 'sfGFP', 'AiiA', 'Aqs1','cellD']
    y0 = [20.0,5.0,5.0, 1.0,1.0, 0.1,0.1,0.2]
    for ara in araList:
        mod_compCir.makeModel(y0,ara,D6,D12,Dc)
        mod_compCir.simulateModel()
        mod_compCir.analysis_plots()
        print(f'Completed ara dose = {ara}')

    # mod_compCir.timeSerPlots()
    # mod_compCir.saveDoseResponse()
    # mod_compCir.plotDoseResponse()

    #end code
    t2= time.time()
    totalSec= t2-t1
    Sec=int(totalSec%60)
    Hrs=int(totalSec//3600)
    Min=int((totalSec%3600)//60)

    print ("Program completed in %sHr:%sMin:%ssec\n"%(Hrs,Min,Sec))
