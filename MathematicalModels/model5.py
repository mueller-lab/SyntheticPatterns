"""
Amit Landge

Goal -To model the complete synthetic patterning circuit

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

    def getLawnRD(self, secR = 0.5, shift= 0.0, uniF = False):
        """
        a method to define lawn-occupied area of the agar
        reaction-part of the system happens only on the lawn occupied area
        diffusion hapeens throught entire agar area

        secR: the fraction of the domain length that is occupied by cells
        shift: to shift the lawn-center to right (+) or left (-) as % of domain-length
        """
        secD = 1. - secR # agar-area fraction without the cell-lawn
        domeShp = np.zeros(self.nx)

        if shift==0.0:
            self.leftInd = int((0.5*secD+shift)*self.nx)
            self.rightInd =int(((1-0.5*secD)+shift)*self.nx)
            # len80 = np.zeros(self.nx)
            # len80[self.leftInd:self.rightInd] = 1.0
            # # len80[:] = 1.0
            # self.len80 = len80.astype('float')

            #midIndex - center of lawn (middle of leftInd and rightInd)
            self.midInd = int((self.leftInd+self.rightInd)/2)

            # to make initial dome-shaped trigger
            domeShape = np.linspace(0.0, np.pi,int(secR*self.nx))
            # domeShape = np.linspace(np.pi/4, 3*np.pi/4,int(secR*self.nx))
            domeShape= 0.2+np.sin(domeShape) # make a domeShape curve using sine wave

            # add the dome-shpae to the initial domeShp array - make sure domeShape has same length as the section of domeShp
            lendS = len(domeShape)
            domeShp[self.leftInd:lendS+self.leftInd] = domeShape
            self.domeShp = domeShp.astype('float')
        else:
            self.leftInd = max(0,int((0.5*secD+shift)*self.nx))
            self.rightInd =int(((1-0.5*secD)+shift)*self.nx)
            # len80 = np.zeros(self.nx)
            # len80[self.leftInd:self.rightInd] = 1.0
            # # len80[:] = 1.0
            # self.len80 = len80.astype('float')

            #midIndex - center of lawn (middle of leftInd and rightInd)
            if shift>(-0.2):
                self.midInd = int((self.leftInd+self.rightInd)/2)
            else:
                self.midInd = self.rightInd - int(6/self.dx)

            # to make initial dome-shaped trigger
            domeShape = np.linspace(3*np.pi/4,np.pi,self.rightInd - self.leftInd)
            # domeShape = np.linspace(np.pi/4, 3*np.pi/4,int(secR*self.nx))
            domeShape= 0.2+np.sin(domeShape) # make a domeShape curve using sine wave

            # add the dome-shpae to the initial domeShp array - make sure domeShape has same length as the section of domeShp
            lendS = len(domeShape)
            domeShp[self.leftInd:lendS+self.leftInd] = domeShape
            self.domeShp = domeShp.astype('float')

        # to test uniform cellD
        if uniF: # set to False by default
            self.domeShp[self.leftInd:self.rightInd] = 1.0

        print(f'left, right, mid indices are = {self.leftInd, self.rightInd, self.midInd}')

    def makeModel(self,y0,ara,D6,D12,Dc):
        """
        To define the system of equations
        """
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
        # parFile.close()

        # to restrict the reactions (cells) in the middle secR % of the domain
        self.getLawnRD(secR = self.secR, shift= self.shift, uniF= self.uniF)

        # starting state
        self.InVal= y0

        # create CellVariables
        self.vL = [] # list of CellVariables
        for i, vName in enumerate(self.vNames):
            var = CellVariable(name = vName, mesh = self.mesh, value = self.InVal[i]*self.domeShp) #
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
        # cap6SecR = 1.0 + 0.1*(1.0 - self.secR)
        self.AHL6cap = (K_cell- self.vL[7])/(self.secR*K_cell)
        self.AHL12cap = (1.1*K_cell- self.vL[7])/(1.1*self.secR*K_cell)

        self.pLux = self.LuxR_AHL6 + self.LuxR_AHL12
        self.pLas = self.LasR_AHL12

        self.eq00 = self.vL[7]*(k_prot*(self.pLux + self.pLas) - kd_ssrA*self.vL[0]) # LuxI

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

        # self.mid0 = int(0.2*self.nx)
        # self.mid1 = int(0.8*self.nx)

        # a faster method to write data to csv file (using pandas)
        self.outDict = {}
        for i, vName in enumerate(self.vNames):
            self.outDict[vName] = []
            self.outDict[vName].append(np.array(self.vL[i]))

    def simulateModel(self, ):
        """
        use eq.solve to interatively solve the system
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

            #append data to outDict
            for i, vName in enumerate(self.vNames):
                self.outDict[vName].append(np.array(self.vL[i]))

            #update step and tNow
            step+=1
            self.tNow+=self.dt

            #keep track of progress
            if step%100 ==0:
                print(f'completed step = {step}')

        # finally saving all data to .csv files (separate file of each vName)
        colNames = [i*self.dx for i in range(self.nx)]
        for i, vName in enumerate(self.vNames):
            df_v = pd.DataFrame(self.outDict[vName], columns = colNames, index = self.tAx)
            fName = f"{vName}.csv"
            df_v.to_csv(self.araPath/fName)

    def analysis_plots(self, subVal = 4, xticInt = 120):
        # dict to store time series of all variables
        timeSerDict={}

        # plot the results as space-time graphs and save
        for j, vName in enumerate(self.vNames):
            fName = f"{vName}.csv"
            fPath = self.araPath/fName
            data_df = pd.read_csv(str(fPath),header=0, index_col = 0)

            vals = data_df.to_numpy(dtype='float')

            # select the data - rows - timepoints, columns- gridpoints (0.5*nx to (1+secR)*0.5*nx)
            # gridMid = int(0.5*self.nx)
            # gridEdg = int((1+self.secR)*gridMid)
            # gridEdg = int((1+0.65)*gridMid)
            vals_sub = vals[::subVal,self.midInd:self.rightInd]

            vals_sub_left = vals[::subVal,self.leftInd:self.midInd]
            vals_sub_left = np.flip(vals_sub_left, 1)

            # for Big lawns limit vals_sub
            if self.secR==1.0:
                vals_sub = vals_sub[:,0:int(5/self.dx)]


            # assign colormap
            if ('mCherry' in vName) or ('AHL12' in vName):
                cmap_kymo = 'magma'
            else:
                cmap_kymo = 'viridis'

            # set min max for ara dose response simulations
            vmin = 0
            if len(self.araList)>1:
                if vName == 'sfGFP':
                    vmax = max(170., np.max(vals_sub))
                elif vName =='mCherry':
                    vmax = max(40., np.max(vals_sub))
                else:
                    vmax = np.max(vals_sub)
            else:
                vmax = np.max(vals_sub)

            #code to make space time plots
            fig, ax  = plt.subplots(nrows = 1, ncols = 1, figsize= (2.0,2.0))
            pos = ax.imshow(vals_sub.T, origin = 'lower', cmap = cmap_kymo, vmin = vmin, vmax = vmax)

            ax.grid(False)
            ax.set_ylabel("r (mm)")
            ax.set_xlabel("Time (h)")

            # set x and y-ticks
            yAx = np.array(data_df.columns)
            yAx = yAx.astype(float)
            yAx = yAx[0:(self.rightInd-self.midInd)]

            # for Big lawns limit yAx -
            if self.secR==1.0:
                yAx = yAx[0:int(5/self.dx)]

            tAx = np.array(data_df.index)
            tAx = tAx.astype(float) # time in min
            tAx_sub = tAx[::subVal]

            tAx_h = tAx_sub/60.0 # time in hr

            yTic = [int(2*i/yAx[1]) for i in range(1+int(yAx[-1]/2))]
            y_labels = [2*i for i in range(1+int(yAx[-1]/2))]
            ax.set_yticks(yTic)
            ax.set_yticklabels(y_labels)

            xTic = [int(xticInt*i/tAx_h[1]) for i in range(1+int(tAx_h[-1]/xticInt))]
            x_labels = [xticInt*i for i in range(1+int(tAx_h[-1]/xticInt))]
            ax.set_xticks(xTic)
            ax.set_xticklabels(x_labels)

            ax.tick_params(axis='both', which='major', labelsize=6)

            # add colorbar
            plt.colorbar(pos, fraction = 0.04, pad = 0.3, orientation = 'horizontal')

            plt.tight_layout()

            #remove top and right spines
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

            figName_0 = f"{vName}_SpaceTime.png"
            figName_1 = f"{vName}_SpaceTime.pdf"

            # plt.show()
            plt.savefig(self.araPath/figName_0, transparent=True, dpi = 300)
            plt.savefig(self.araPath/figName_1, transparent=True)
            plt.close()

            #######################################################
            if -0.2<self.shift<0.0:
                # plot the profile from the left of lawn center (flipped)
                #code to make space time plots
                fig, ax  = plt.subplots(nrows = 1, ncols = 1, figsize= (2.0,2.0))
                pos = ax.imshow(vals_sub_left.T, origin = 'lower', cmap = cmap_kymo)

                ax.grid(False)
                ax.set_ylabel("r (mm)")
                ax.set_xlabel("Time (h)")

                # set x and y-ticks
                yAx = np.array(data_df.columns)
                yAx = yAx.astype(float)
                yAx = yAx[0:(self.midInd-self.leftInd)]

                yTic = [int(2*i/yAx[1]) for i in range(1+int(yAx[-1]/2))]
                y_labels = [2*i for i in range(1+int(yAx[-1]/2))]
                ax.set_yticks(yTic)
                ax.set_yticklabels(y_labels)

                ax.set_xticks(xTic)
                ax.set_xticklabels(x_labels)

                ax.tick_params(axis='both', which='major', labelsize=6)

                # add colorbar
                plt.colorbar(pos, fraction = 0.04, pad = 0.3, orientation = 'horizontal')

                plt.tight_layout()

                #remove top and right spines
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)

                figName_0 = f"{vName}_SpaceTime_left.png"
                figName_1 = f"{vName}_SpaceTime_left.pdf"

                # plt.show()
                plt.savefig(self.araPath/figName_0, transparent=True, dpi = 300)
                plt.savefig(self.araPath/figName_1, transparent=True)
                plt.close()

            # get an array of mean values from the cell lawn area(10% of the total region)
            # secR = 0.65 # to only get values equivalent to 'Original'
            # secD = 1. - secR # agar-area fraction without the cell-lawn
            # self.leftInd = max(0, int((0.5*secD+shift)*self.nx))
            # self.rightInd = min(self.nx, int(((1-0.5*secD)+shift)*self.nx))

            meanCent = np.mean(vals[:,self.leftInd:self.rightInd], axis = 1)

            # add timeseries to the dict
            timeSerDict[vName]= meanCent

            # add final mCherry conc. to a dict to make the dose response curve
            if vName =="mCherry":
                if len(meanCent)>int(360/self.dt):
                    self.doseRes[ara] = meanCent[int(360/self.dt)]
                else:
                    self.doseRes[ara] = meanCent[-1]

        # save time series of all variable in the respective ahl dir
        timeSerDF = pd.DataFrame.from_dict(timeSerDict)
        timeSerDF.to_csv(self.araPath/'timeSerAll.csv')

    def timeSerPlots(self,):
        self.tAxH = np.array(self.tAx)/60
        for vPlot in self.vNames:
            # magma for mCherry
            if 'mCherry' in vPlot:
                # make a color dict to assign colors to typeId
                viridis = matplotlib.colormaps.get_cmap('magma')
                cmap0 = viridis(np.linspace(0, 1, 1+len(self.araList)))
            else:
                # make a color dict to assign colors to typeId
                viridis = matplotlib.colormaps.get_cmap('viridis')
                cmap0 = viridis(np.linspace(0, 1, len(self.araList)))

            # load var timeseries at all ahl treatments
            normSeries = []
            for ara in self.araList:
                fName = self.resultPath/f"ara_{ara:.2f}"/'timeSerAll.csv'
                data_df = pd.read_csv(str(fName),header=0, index_col = 0)
                vals = np.array(data_df[vPlot])
                normSeries.append(vals)

            normArr = np.array(normSeries)
            # to make time-series plots with mean values from domain center
            fig, ax  = plt.subplots(nrows = 1, ncols = 1, figsize= (2.4,1.6))

            for i, ara in enumerate(self.araList):
                ax.plot(self.tAxH, normArr[i,:], linewidth = 1, c=cmap0[i], label = f'{ara:.2e}')

            ax.set_xlabel("Time (h)")
            ax.set_ylabel("Output (A.U.)")
            ax.set_ylim(0.0,1.02*np.amax(normArr))
            ax.set_xlim(0.0,1.02*self.tAxH[-1])
            ax.tick_params(axis='both', which='major', labelsize=6)
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
            ax.get_legend().set_title("Arabinose (nM)")

            #remove top and right spines
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            plt.tight_layout()

            figName_0 = f"TimeSeries_{vPlot}.png"
            figName_1 = f"TimeSeries_{vPlot}.pdf"

            # plt.show()
            plt.savefig(self.resultPath/figName_0, transparent=True, dpi = 300)
            plt.savefig(self.resultPath/figName_1, transparent=True)
            plt.close()

            # to make time-series plots with mean values from domain center
            normArr0 = normArr/np.amax(normArr)
            fig, ax  = plt.subplots(nrows = 1, ncols = 1, figsize= (2.4,1.6))

            for i, ara in enumerate(self.araList):
                ax.plot(self.tAxH, normArr0[i,:], linewidth = 1, c=cmap0[i], label = f'{ara:.2e}')

            ax.set_xlabel("Time (h)")
            ax.set_ylabel("Output (A.U.)")
            ax.set_ylim(0.0,1.02)
            ax.set_xlim(0.0,1.02*self.tAxH[-1])
            ax.tick_params(axis='both', which='major', labelsize=6)
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
            ax.get_legend().set_title("Arabinose (nM)")

            #remove top and right spines
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            plt.tight_layout()

            figName_0 = f"normTimeSeries_{vPlot}.png"
            figName_1 = f"normTimeSeries_{vPlot}.pdf"

            # plt.show()
            plt.savefig(self.resultPath/figName_0, transparent=True, dpi = 300)
            plt.savefig(self.resultPath/figName_1, transparent=True)
            plt.close()

    def saveDoseResponse(self,):
        # make data array
        doseRes = np.array(list(self.doseRes.items()))

        #save doseRes dict as a pandas DataFrame
        doseDf = pd.DataFrame(self.doseRes.items(), columns = ['ara', 'mCherry'])
        doseDf.to_csv(self.resultPath/"doseRes.csv")

    def plotDoseResponse(self,):
        # to make dose-response plots with mean values from domain center
        fig, ax  = plt.subplots(nrows = 1, ncols = 1, figsize= (1.6,1.6))

        # read data
        dataDF = pd.read_csv(self.resultPath/"doseRes.csv", header=0, index_col = 0)
        doseRes = np.array(dataDF)
        xAra = doseRes[:,0]
        yRes = doseRes[:,1]

        # normalize mCherry output to 1
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

    def plotMeanDyn(self, ):
        """
        to plot mCherry and sfGFP over time similar to Fig. 4e
        """
        # # make colormap
        viridis = matplotlib.colormaps.get_cmap('viridis')
        cmapV = viridis(np.linspace(0, 1, 5))

        magma = matplotlib.colormaps.get_cmap('magma')
        cmapM = magma(np.linspace(0, 1, 5))

        # make a dict to assign cmap to ch
        clrDict = {'sfGFP':cmapV[3], 'mCherry':cmapM[3]}
        for ara in self.araList:
            self.araPath=self.resultPath/f"ara_{ara:.2f}"
            # make figure
            fig, ax  = plt.subplots(nrows = 1, ncols = 1, figsize= (2.4,1.6))
            for vPlot in ['sfGFP', 'mCherry']:
                fName = self.araPath/'timeSerAll.csv'
                data_df = pd.read_csv(str(fName),header=0, index_col = 0)
                vals = np.array(data_df[vPlot])
                vals = vals/np.max(vals) # scale to 1
                ax.plot(self.tAxH, vals, marker= ".", markersize=1, color =clrDict[vPlot], label = vPlot, linewidth = 1)

            ax.grid(False)
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

            ax.set_xlabel("Time (h)")
            ax.set_ylabel("Mean intensity (A.U.)")
            # ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

            ax.set_ylim([0, 1.02])
            ax.set_xlim(self.tAxH[0], self.tAxH[-1])

            ax.tick_params(axis='both', which='major', labelsize=6)

            plt.tight_layout()

            #remove top and right spines
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

            figName_0 = f"mean_vs_time.png"
            figName_1 = f"mean_vs_time.pdf"

            plt.savefig(self.araPath/figName_0, transparent=True, dpi = 300)
            plt.savefig(self.araPath/figName_1, transparent=True)
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
    plt.rcParams['pdf.fonttype'] = 42

    #get current dir and datetime
    cwdPath=Path(os.path.abspath(os.getcwd()))
    now=datetime.now()
    datetime_str=(now.strftime("%Y%m%d_")) # %H%M%S_

    #make output directory
    dirName=str(datetime_str+ "mod5_compCir_Big") # the positive loop model

    resultPath=cwdPath/'data_1D'/dirName
    resultPath.mkdir(mode=0o777, parents=True, exist_ok=True)
    print(f"Created result directory {dirName} at {time.time() - t1} sec ...")

    #define global parameters
    steps = 360
    dt= 4.0 # time in min
    gridSize = 16 # mm
    dx= 0.05 # mm
    nx = int(gridSize/dx)
    vNames = ['LuxI', 'AHL12', 'AHL6', 'mCherry', 'sfGFP', 'AiiA', \
    'Aqs1','cellD']

    # define diffusion rates (per min)
    D6= 0.0 # 2.1*1E-2 # AHL6 diffusion rate(act)
    D12 = 0.0 # 2.1*1E-2 # AHL12 diffusion rate (inh)

    # synthesis rates (1/min, )
    k_prot = 0.2 # synthesis time  ~ 5 min
    k_gfp = 0.1 # synthesis + maturation time ~ 10 min
    k_mChr = 0.02 # synthesis + maturation time ~ 50 min
    k_ahl = 0.002 # AHL synthesis rate (by LuxI)

    # degradation rates (unit = 1/min, )
    # note: half-life of 70 min --> degradation rate = 0.01 1/min
    kd_ssrA = 0.018 # half life of 38.5 min - based on http://parts.igem.org/Part:BBa_K1399004
    kd_prot = 0.0056 #t1/2 ~ 2 h
    kd_ahl = 0.0011 # ahl decay (hydroslysis, loss into agar), half-life ~ 10 h
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

    secR = 1.0 #0.65 # area-fraction of agar occupied by cell lawn (small - 0.5, big - 1.0)
    shift = 0.0 # a parameter to shift the lawn center to left, as fraction of total length, value between (-0.9 to 0)
    uniF = False # True

    # test uniform-conditions (diffusion = 0 )
    # D6, D12, Dc = [0.,0.,0.]

    # arabinose treatments
    araList = [0.0] # ,  1., 1E3, 1E6

    mod_compCir = RD_simulate(nx, dx, dt, steps, vNames, resultPath, araList, secR, shift, uniF)

    # ['LuxI', 'AHL12', 'AHL6', 'mCherry', 'sfGFP', 'AiiA', 'Aqs1','cellD']
    y0 = [20.0,5.0,5.0, 1.0,1.0, 0.1,0.1,0.2]
    for ara in araList:
        mod_compCir.makeModel(y0,ara,D6,D12,Dc)
        mod_compCir.simulateModel()
        mod_compCir.analysis_plots(subVal = 6, xticInt = 10)
        print(f'Completed ara dose = {ara}')

    mod_compCir.timeSerPlots()
    mod_compCir.plotMeanDyn()

    mod_compCir.saveDoseResponse()
    # mod_compCir.plotDoseResponse()

    #end code
    t2= time.time()
    totalSec= t2-t1
    Sec=int(totalSec%60)
    Hrs=int(totalSec//3600)
    Min=int((totalSec%3600)//60)

    print ("Program completed in %sHr:%sMin:%ssec\n"%(Hrs,Min,Sec))
