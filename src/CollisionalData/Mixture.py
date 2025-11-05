#!/usr/bin/env python3
import sys
from mpi4py import MPI
sys.path.append(sys.path[0] + '/../../build')
import muffin				        ## module in C++ doing the simulation
import glob				            ## module for unix pathname pattern expansioni
import matplotlib.pyplot as plt 	## plots 
import os 				            ## operating system interface
import numpy as np			        ## Mathematic operations
import shutil                       ## Remove folder
import scipy.constants as phy_const ## Physical constants
from scipy.optimize import curve_fit
sys.path.append('./')
import CrossSection as CS

class Mixture:

    def __init__(self, name=None, electronData=None, ionData=None):
        self.name         = name
        self.electronData = electronData
        self.ionData      = ionData
        self.Temperature  = np.linspace(1.,20,1000)
        self.Delta        = []
        self.Ktotal       = np.linspace(1.5,20,1000)
        self.Losses       = np.zeros(1000)
        self.Omega11_en   = np.zeros(1000)
        self.Omega12_en   = np.zeros(1000)
        self.Omega13_en   = np.zeros(1000)
        self.Omega14_en   = np.zeros(1000)
        self.Omega22_ee   = np.zeros(1000)
        self.Omega23_ee   = np.zeros(1000)
        self.Omega24_ee   = np.zeros(1000)
        self.CoeffsLosses = []
        self.CoeffsLossesDelta = []
        self.CoeffsOmega11_en = []
        self.CoeffsOmega12_en = []
        self.CoeffsOmega13_en = []
        self.CoeffsOmega14_en = []
        self.CoeffsOmega22_ee = []
        self.CoeffsOmega23_ee = []
        self.CoeffsOmega24_ee = []
        
    def computeFullLosses(self):
        self.Losses      = np.zeros(1000)
        for collision in self.electronData:
            collision.computeRate()
            self.Losses[:] += collision.potential*collision.K[:]
            
    def computeFullLossesDelta(self):
        self.LossesDelta    = np.zeros((101, 1000))
        for collision in self.electronData:
            collision.computeRateDeltaZero()
            self.Delta = collision.Delta
            for iDelta, Delta in enumerate(collision.Delta):
                for iTemp, Temp in enumerate(collision.Temperature):
                    self.LossesDelta[iDelta,iTemp] += collision.KDeltaZero[iDelta,iTemp]*collision.potential
        
    def computeFullLossesFourthDelta(self):
        self.LossesFourthDelta    = np.zeros((101, 1000))
        for collision in self.electronData:
            collision.computeRateDeltaZero()
            collision.computeRateDeltaOne()
            self.Delta = collision.Delta
            for iDelta, Delta in enumerate(collision.Delta):
                for iTemp, Temp in enumerate(collision.Temperature):
                    self.LossesFourthDelta[iDelta, iTemp] += collision.KDeltaZero[iDelta,iTemp]*collision.potential**2/collision.Temperature[iTemp]**2 +  0.5*collision.KDeltaOne[iDelta,iTemp]*collision.potential/collision.Temperature[iTemp]
                    #               collision.KDeltaOne[iDelta,iTemp]*collision.potential/collision.Temperature[iTemp]
                    # collision.KDeltaZero[iDelta,iTemp]*collision.potential**2/collision.Temperature[iTemp]**2 +
            
    def computeOmega11_en(self):
        self.Omega11_en      = np.zeros(1000)
        for collision in self.electronData:
            if(collision.mechanism == "e-neutral"):
                collision.computeOmega11()
                self.Omega11_en[:] += collision.Omega11[:]
    
    def computeOmega12_en(self):
        self.Omega12_en      = np.zeros(1000)
        for collision in self.electronData:
            if(collision.mechanism == "e-neutral"):
                collision.computeOmega12()
                self.Omega12_en[:] += collision.Omega12[:]
    
    def computeOmega13_en(self):
        self.Omega13_en      = np.zeros(1000)
        for collision in self.electronData:
            if(collision.mechanism == "e-neutral"):
                collision.computeOmega13()
                self.Omega13_en[:] += collision.Omega13[:]
    
    def computeOmega14_en(self):
        self.Omega14_en      = np.zeros(1000)
        for collision in self.electronData:
            if(collision.mechanism == "e-neutral"):
                collision.computeOmega14()
                self.Omega14_en[:] += collision.Omega14[:]
    
    def computeOmega22_ee(self):
        self.Omega22_ee      = np.zeros(1000)
        for collision in self.electronData:
            if(collision.mechanism == "e-e"):
                collision.computeOmega22()
                self.Omega22_ee[:] = collision.Omega22[:]
    
    def computeOmega23_ee(self):
        self.Omega23_ee      = np.zeros(1000)
        for collision in self.electronData:
            if(collision.mechanism == "e-e"):
                collision.computeOmega23()
                self.Omega23_ee[:] = collision.Omega23[:]
    
    def computeOmega24_ee(self):
        self.Omega24_ee      = np.zeros(1000)
        for collision in self.electronData:
            if(collision.mechanism == "e-e"):
                collision.computeOmega24()
                self.Omega24_ee[:] = collision.Omega24[:]
        
    def fitFullLosses(self):
        Approx            = curve_fit(lambda t,a,b,c,d,f: a*t**f*np.exp(-b*t + c*t**2 - d*t**3.),  1./self.Temperature,  self.Losses,  p0=(self.Losses[-1], 0.,0.,0.,0.))
        self.CoeffsLosses = curve_fit(lambda t,a,b,c,d,e,f: a*t**f*np.exp(-b*t+ c*t**2 - d*t**3. + e*t**4),  1./self.Temperature,  self.Losses,  p0=(Approx[0][0], Approx[0][1], Approx[0][2], Approx[0][3], 0., Approx[0][-1]), method='dogbox')[0]
        
#    def fitFullLossesDelta(self):
#        for iDelta in range(0,101):
#        ## TODO
#            print(iDelta)
#            plt.plot(self.Temperature, np.log(self.LossesDelta[iDelta,:]))
#        plt.show()
#
#            Approx            = curve_fit(lambda t,a,b,c,d,f: a*t**f*np.exp(-b*t + c*t**2 - d*t**3.),  1./self.Temperature,  self.LossesDelta[iDelta,:],  p0=(self.LossesDelta[iDelta,-1], 0.,0.,0.,0.))
#            CoeffsLosses_Delta = curve_fit(lambda t,a,b,c,d,e,f: a*t**f*np.exp(-b*t+ c*t**2 - d*t**3. + e*t**4),  1./self.Temperature,  self.LossesDelta[iDelta,:],  p0=(Approx[0][0], Approx[0][1], Approx[0][2], Approx[0][3], 0., Approx[0][-1]), method='dogbox')[0]
#            self.CoeffsLossesDelta.append(CoeffsLosses_Delta)
#        print(self.CoeffsLossesDelta)
    
    def fitOmega11_en(self):
        Approx            = curve_fit(lambda t,a,b,c,d,f: a*t**f*np.exp(-b*t + c*t**2 - d*t**3.),  1./self.Temperature,  self.Omega11_en,  p0=(self.Omega11_en[-1], 0.,0.,0.,0.))
        self.CoeffsOmega11_en = curve_fit(lambda t,a,b,c,d,e,f: a*t**f*np.exp(-b*t+ c*t**2 - d*t**3. + e*t**4),  1./self.Temperature,  self.Omega11_en,  p0=(Approx[0][0], Approx[0][1], Approx[0][2], Approx[0][3], 0., Approx[0][-1]), method='dogbox')[0]
    
    def fitOmega12_en(self):
        Approx            = curve_fit(lambda t,a,b,c,d,f: a*t**f*np.exp(-b*t + c*t**2 - d*t**3.),  1./self.Temperature,  self.Omega12_en,  p0=(self.Omega12_en[-1], 0.,0.,0.,0.))
        self.CoeffsOmega12_en = curve_fit(lambda t,a,b,c,d,e,f: a*t**f*np.exp(-b*t+ c*t**2 - d*t**3. + e*t**4),  1./self.Temperature,  self.Omega12_en,  p0=(Approx[0][0], Approx[0][1], Approx[0][2], Approx[0][3], 0., Approx[0][-1]), method='dogbox')[0]
        
    def fitOmega13_en(self):
        Approx            = curve_fit(lambda t,a,b,c,d,f: a*t**f*np.exp(-b*t + c*t**2 - d*t**3.),  1./self.Temperature,  self.Omega13_en,  p0=(self.Omega13_en[-1], 0.,0.,0.,0.))
        self.CoeffsOmega13_en = curve_fit(lambda t,a,b,c,d,e,f: a*t**f*np.exp(-b*t+ c*t**2 - d*t**3. + e*t**4),  1./self.Temperature,  self.Omega13_en,  p0=(Approx[0][0], Approx[0][1], Approx[0][2], Approx[0][3], 0., Approx[0][-1]), method='dogbox')[0]
        
    def fitOmega14_en(self):
        Approx            = curve_fit(lambda t,a,b,c,d,f: a*t**f*np.exp(-b*t + c*t**2 - d*t**3.),  1./self.Temperature,  self.Omega14_en,  p0=(self.Omega14_en[-1], 0.,0.,0.,0.))
        self.CoeffsOmega14_en = curve_fit(lambda t,a,b,c,d,e,f: a*t**f*np.exp(-b*t+ c*t**2 - d*t**3. + e*t**4),  1./self.Temperature,  self.Omega14_en,  p0=(Approx[0][0], Approx[0][1], Approx[0][2], Approx[0][3], 0., Approx[0][-1]), method='dogbox')[0]
        
    def fitOmega22_ee(self):
        Approx            = curve_fit(lambda t,a,b,c,d,f: a*t**f*np.exp(-b*t + c*t**2 - d*t**3.),  1./self.Temperature,  self.Omega22_ee,  p0=(self.Omega22_ee[-1], 0.,0.,0.,0.))
        self.CoeffsOmega22_ee = curve_fit(lambda t,a,b,c,d,e,f: a*t**f*np.exp(-b*t+ c*t**2 - d*t**3. + e*t**4),  1./self.Temperature,  self.Omega22_ee,  p0=(Approx[0][0], Approx[0][1], Approx[0][2], Approx[0][3], 0., Approx[0][-1]), method='dogbox')[0]
        
    def fitOmega23_ee(self):
        Approx            = curve_fit(lambda t,a,b,c,d,f: a*t**f*np.exp(-b*t + c*t**2 - d*t**3.),  1./self.Temperature,  self.Omega23_ee,  p0=(self.Omega23_ee[-1], 0.,0.,0.,0.))
        self.CoeffsOmega23_ee = curve_fit(lambda t,a,b,c,d,e,f: a*t**f*np.exp(-b*t+ c*t**2 - d*t**3. + e*t**4),  1./self.Temperature,  self.Omega23_ee,  p0=(Approx[0][0], Approx[0][1], Approx[0][2], Approx[0][3], 0., Approx[0][-1]), method='dogbox')[0]
        
    def fitOmega24_ee(self):
        Approx            = curve_fit(lambda t,a,b,c,d,f: a*t**f*np.exp(-b*t + c*t**2 - d*t**3.),  1./self.Temperature,  self.Omega24_ee,  p0=(self.Omega24_ee[-1], 0.,0.,0.,0.))
        self.CoeffsOmega24_ee = curve_fit(lambda t,a,b,c,d,e,f: a*t**f*np.exp(-b*t+ c*t**2 - d*t**3. + e*t**4),  1./self.Temperature,  self.Omega24_ee,  p0=(Approx[0][0], Approx[0][1], Approx[0][2], Approx[0][3], 0., Approx[0][-1]), method='dogbox')[0]
    
    def plotFitFullLosses(self, save=True):
        def fit_function(T):
            a = self.CoeffsLosses[0]
            b = self.CoeffsLosses[1]
            c = self.CoeffsLosses[2]
            d = self.CoeffsLosses[3]
            e = self.CoeffsLosses[4]
            f = self.CoeffsLosses[5]

            t = (1/T)

            return a*t**f*np.exp(-b*t+ c*t**2 - d*t**3. + e*t**4)
            
        f, ax = plt.subplots()
        fit = fit_function(self.Temperature)
        ax.plot(self.Temperature, fit, label = 'fit', linestyle='None', marker='x', markevery=15, color='k')
        ax.plot(self.Temperature, self.Losses, label = 'elec. losses', color='k')
        
        ax.set_ylabel(r'$K^*\cdot \phi^*$ $[eV\cdot m^{3}/s]$')
        ax.set_xlabel('T [eV]')
        ax.set_xlim([2,20])
        ax.set_yscale('log')

        plt.legend(loc = 'best', fontsize = 5)
        plt.title('Electron losses fit')
        
        ax.grid(True)
        file_name = "./ElecLosses"
        if save:
            print("Saving the figure of the electron losses in ./ElecLosses.pdf")
            #if path.exists(file_name+".pdf"):
            #    file_name = file_name+"_1"
            file_name = file_name+".pdf"
            plt.savefig(file_name, bbox_inches='tight')
            plt.close()
        else:
            plt.show()
            plt.close()
            
    def plotFitOmega11_en(self, save=True):
        def fit_function(T):
            a = self.CoeffsOmega11_en[0]
            b = self.CoeffsOmega11_en[1]
            c = self.CoeffsOmega11_en[2]
            d = self.CoeffsOmega11_en[3]
            e = self.CoeffsOmega11_en[4]
            f = self.CoeffsOmega11_en[5]

            t = (1/T)

            return a*t**f*np.exp(-b*t+ c*t**2 - d*t**3. + e*t**4)
            
        f, ax = plt.subplots()
        fit = fit_function(self.Temperature)
        ax.plot(self.Temperature, fit, label = 'fit', linestyle='None', marker='x', markevery=15, color='k')
        ax.plot(self.Temperature, self.Omega11_en, label = r'total $\Omega^{(1,1)}_{en}$', color='k')
        
        ax.set_ylabel(r'$\Omega^{(1,1)}_{en}$ $[m^{3}/s]$')
        ax.set_xlabel('T [eV]')
        ax.set_xlim([2,20])
        ax.set_yscale('log')

        plt.legend(loc = 'best', fontsize = 5)
        plt.title('Total $\Omega^{(1,1)}_{en}$ fit')
        
        ax.grid(True)
        file_name = "./Omega11_en"
        if save:
            print("Saving the figure of the Omega11_en in ./Omega11_en.pdf")
            #if path.exists(file_name+".pdf"):
            #    file_name = file_name+"_1"
            file_name = file_name+".pdf"
            plt.savefig(file_name, bbox_inches='tight')
            plt.close()
        else:
            plt.show()
            plt.close()
    
    def plotFitOmega12_en(self, save=True):
        def fit_function(T):
            a = self.CoeffsOmega12_en[0]
            b = self.CoeffsOmega12_en[1]
            c = self.CoeffsOmega12_en[2]
            d = self.CoeffsOmega12_en[3]
            e = self.CoeffsOmega12_en[4]
            f = self.CoeffsOmega12_en[5]

            t = (1/T)

            return a*t**f*np.exp(-b*t+ c*t**2 - d*t**3. + e*t**4)
            
        f, ax = plt.subplots()
        fit = fit_function(self.Temperature)
        ax.plot(self.Temperature, fit, label = 'fit', linestyle='None', marker='x', markevery=15, color='k')
        ax.plot(self.Temperature, self.Omega12_en, label = r'total $\Omega^{(1,2)}_{en}$', color='k')
        
        ax.set_ylabel(r'$\Omega^{(1,1)}_{en}$ $[m^{3}/s]$')
        ax.set_xlabel('T [eV]')
        ax.set_xlim([2,20])
        ax.set_yscale('log')

        plt.legend(loc = 'best', fontsize = 5)
        plt.title('Total $\Omega^{(1,2)}_{en}$ fit')
        
        ax.grid(True)
        file_name = "./Omega12_en"
        if save:
            print("Saving the figure of the Omega12_en in ./Omega12_en.pdf")
            #if path.exists(file_name+".pdf"):
            #    file_name = file_name+"_1"
            file_name = file_name+".pdf"
            plt.savefig(file_name, bbox_inches='tight')
            plt.close()
        else:
            plt.show()
            plt.close()
            
    def plotFitOmega13_en(self, save=True):
        def fit_function(T):
            a = self.CoeffsOmega13_en[0]
            b = self.CoeffsOmega13_en[1]
            c = self.CoeffsOmega13_en[2]
            d = self.CoeffsOmega13_en[3]
            e = self.CoeffsOmega13_en[4]
            f = self.CoeffsOmega13_en[5]

            t = (1/T)

            return a*t**f*np.exp(-b*t+ c*t**2 - d*t**3. + e*t**4)
            
        f, ax = plt.subplots()
        fit = fit_function(self.Temperature)
        ax.plot(self.Temperature, fit, label = 'fit', linestyle='None', marker='x', markevery=15, color='k')
        ax.plot(self.Temperature, self.Omega13_en, label = r'total $\Omega^{(1,1)}_{en}$', color='k')
        
        ax.set_ylabel(r'$\Omega^{(1,3)}_{en}$ $[m^{3}/s]$')
        ax.set_xlabel('T [eV]')
        ax.set_xlim([2,20])
        ax.set_yscale('log')

        plt.legend(loc = 'best', fontsize = 5)
        plt.title('Total $\Omega^{(1,3)}_{en}$ fit')
        
        ax.grid(True)
        file_name = "./Omega13_en"
        if save:
            print("Saving the figure of the Omega13_en in ./Omega13_en.pdf")
            #if path.exists(file_name+".pdf"):
            #    file_name = file_name+"_1"
            file_name = file_name+".pdf"
            plt.savefig(file_name, bbox_inches='tight')
            plt.close()
        else:
            plt.show()
            plt.close()
            
    def plotFitOmega14_en(self, save=True):
        def fit_function(T):
            a = self.CoeffsOmega14_en[0]
            b = self.CoeffsOmega14_en[1]
            c = self.CoeffsOmega14_en[2]
            d = self.CoeffsOmega14_en[3]
            e = self.CoeffsOmega14_en[4]
            f = self.CoeffsOmega14_en[5]

            t = (1/T)

            return a*t**f*np.exp(-b*t+ c*t**2 - d*t**3. + e*t**4)
            
        f, ax = plt.subplots()
        fit = fit_function(self.Temperature)
        ax.plot(self.Temperature, fit, label = 'fit', linestyle='None', marker='x', markevery=15, color='k')
        ax.plot(self.Temperature, self.Omega14_en, label = r'total $\Omega^{(1,4)}_{en}$', color='k')
        
        ax.set_ylabel(r'$\Omega^{(1,4)}_{en}$ $[m^{3}/s]$')
        ax.set_xlabel('T [eV]')
        ax.set_xlim([2,20])
        ax.set_yscale('log')

        plt.legend(loc = 'best', fontsize = 5)
        plt.title('Total $\Omega^{(1,4)}_{en}$ fit')
        
        ax.grid(True)
        file_name = "./Omega14_en"
        if save:
            print("Saving the figure of the Omega14_en in ./Omega14_en.pdf")
            #if path.exists(file_name+".pdf"):
            #    file_name = file_name+"_1"
            file_name = file_name+".pdf"
            plt.savefig(file_name, bbox_inches='tight')
            plt.close()
        else:
            plt.show()
            plt.close()
    
    
    def plotFitOmega22_ee(self, save=True):
        def fit_function(T):
            a = self.CoeffsOmega22_ee[0]
            b = self.CoeffsOmega22_ee[1]
            c = self.CoeffsOmega22_ee[2]
            d = self.CoeffsOmega22_ee[3]
            e = self.CoeffsOmega22_ee[4]
            f = self.CoeffsOmega22_ee[5]

            t = (1/T)

            return a*t**f*np.exp(-b*t+ c*t**2 - d*t**3. + e*t**4)
            
        f, ax = plt.subplots()
        fit = fit_function(self.Temperature)
        ax.plot(self.Temperature, fit, label = 'fit', linestyle='None', marker='x', markevery=15, color='k')
        ax.plot(self.Temperature, self.Omega22_ee, label = r'total $\Omega^{(221)}_{ee}$', color='k')
        
        ax.set_ylabel(r'$\Omega^{(2,2)}_{ee}$ $[m^{3}/s]$')
        ax.set_xlabel('T [eV]')
        ax.set_xlim([2,20])
        ax.set_yscale('log')

        plt.legend(loc = 'best', fontsize = 5)
        plt.title('Total $\Omega^{(2,2)}_{ee}$ fit')
        
        ax.grid(True)
        file_name = "./Omega22_ee"
        if save:
            print("Saving the figure of the Omega22_ee in ./Omega22_ee.pdf")
            file_name = file_name+".pdf"
            plt.savefig(file_name, bbox_inches='tight')
            plt.close()
        else:
            plt.show()
            plt.close()
    
    def plotFitOmega23_ee(self, save=True):
        def fit_function(T):
            a = self.CoeffsOmega23_ee[0]
            b = self.CoeffsOmega23_ee[1]
            c = self.CoeffsOmega23_ee[2]
            d = self.CoeffsOmega23_ee[3]
            e = self.CoeffsOmega23_ee[4]
            f = self.CoeffsOmega23_ee[5]

            t = (1/T)

            return a*t**f*np.exp(-b*t+ c*t**2 - d*t**3. + e*t**4)
            
        f, ax = plt.subplots()
        fit = fit_function(self.Temperature)
        ax.plot(self.Temperature, fit, label = 'fit', linestyle='None', marker='x', markevery=15, color='k')
        ax.plot(self.Temperature, self.Omega23_ee, label = r'total $\Omega^{(2,3)}_{ee}$', color='k')
        
        ax.set_ylabel(r'$\Omega^{(2,3)}_{ee}$ $[m^{3}/s]$')
        ax.set_xlabel('T [eV]')
        ax.set_xlim([2,20])
        ax.set_yscale('log')

        plt.legend(loc = 'best', fontsize = 5)
        plt.title('Total $\Omega^{(2,3)}_{ee}$ fit')
        
        ax.grid(True)
        file_name = "./Omega23_ee"
        if save:
            print("Saving the figure of the Omega23_ee in ./Omega23_ee.pdf")
            file_name = file_name+".pdf"
            plt.savefig(file_name, bbox_inches='tight')
            plt.close()
        else:
            plt.show()
            plt.close()
            
    def plotFitOmega24_ee(self, save=True):
        def fit_function(T):
            a = self.CoeffsOmega24_ee[0]
            b = self.CoeffsOmega24_ee[1]
            c = self.CoeffsOmega24_ee[2]
            d = self.CoeffsOmega24_ee[3]
            e = self.CoeffsOmega24_ee[4]
            f = self.CoeffsOmega24_ee[5]

            t = (1/T)

            return a*t**f*np.exp(-b*t+ c*t**2 - d*t**3. + e*t**4)
            
        f, ax = plt.subplots()
        fit = fit_function(self.Temperature)
        ax.plot(self.Temperature, fit, label = 'fit', linestyle='None', marker='x', markevery=15, color='k')
        ax.plot(self.Temperature, self.Omega24_ee, label = r'total $\Omega^{(2,4)}_{ee}$', color='k')
        
        ax.set_ylabel(r'$\Omega^{(2,3)}_{ee}$ $[m^{3}/s]$')
        ax.set_xlabel('T [eV]')
        ax.set_xlim([2,20])
        ax.set_yscale('log')

        plt.legend(loc = 'best', fontsize = 5)
        plt.title('Total $\Omega^{(2,4)}_{ee}$ fit')
        
        ax.grid(True)
        file_name = "./Omega24_ee"
        if save:
            print("Saving the figure of the Omega24_ee in ./Omega24_ee.pdf")
            file_name = file_name+".pdf"
            plt.savefig(file_name, bbox_inches='tight')
            plt.close()
        else:
            plt.show()
            plt.close()

    def plotCrossSections(self, save = True):
        import os.path
        from os import path
        # Plot the electron data
        f, ax = plt.subplots()
        
        for collision in self.electronData:
            collision.plotSigma(ax)
            
        ax.set_ylabel(r'$\sigma$ $[m^{2}]$')
        ax.set_xlabel('E [eV]')
        ax.set_yscale('log')
        ax.set_xscale('log')
        #ax.set_ylim([1e-25, 1e-18])
        #ax.set_xlim([1e-2, 1000])
        plt.legend(loc = 'best', fontsize = 5)
        plt.title('Electron collisional Cross-section')
        
        ax.grid(True)
        file_name = "./CrossSections_elec"
        if save:
            print("Saving the figure of the electron cross sections in ./CrossSections_*.pdf")
            #if path.exists(file_name+".pdf"):
            #    file_name = file_name+"_1"
            file_name = file_name+".pdf"
            plt.savefig(file_name, bbox_inches='tight')
            plt.close()
        else:
            plt.show()
            plt.close()
        
        # Plot the ion data
        f, ax = plt.subplots()
        for collision in self.ionData:
            collision.plotSigma(ax)
        
        ax.set_ylabel(r'$\sigma$ $[m^{2}]$')
        ax.set_xlabel('E [eV]')
        ax.set_yscale('log')
        ax.set_xscale('log')
        #ax.set_ylim([1e-25, 1e-18])
        #ax.set_xlim([1e-2, 1000])
        plt.legend(loc = 'best', fontsize = 10)
        plt.title('Ion collisional Cross-section')
        
        ax.grid(True)
        file_name = "./CrossSections_ion"
        if save:
            print("Saving the figure of the ion cross sections in ./CrossSections_*.pdf")
            #if path.exists(file_name+".pdf"):
            #    file_name = file_name+"_1"
            file_name = file_name+".pdf"
            plt.savefig(file_name, bbox_inches='tight')
            plt.close()
        else:
            plt.show()
            plt.close()
        
    def plotRateK(self, save = True):
        print("Saving the figure of the electron Rates in ./Rate_elec.pdf")
        # Plot the data
        f, ax = plt.subplots()

        for collision in self.electronData:
            collision.computeRate()
            collision.plotRate(ax)
            
        ax.set_ylabel(r'$K$ $[m^{3}/s]$')
        ax.set_xlabel(r'$T_e$ [eV]')
        ax.set_yscale('log')
        ax.set_xlim([2., 10])
        #ax.set_ylim([1e-16, 1e-12])
        #plt.legend(loc = 'best', fontsize = 10)
        ax.grid(True)
        
        if save:
            plt.savefig("./Rate_elec.pdf", bbox_inches='tight')
        else:
            plt.show()

    def plotFitRateK(self, save = True):
        # Plot the data
        f, ax = plt.subplots()

        for collision in self.electronData:
            collision.computeRate()
            collision.plotRate(ax)
            collision.fitRate()
            collision.plotFitRate(ax)
            
        ax.set_ylabel(r'$K$ $[m^{3}/s]$')
        ax.set_xlabel(r'$T_e$ [eV]')
        ax.set_yscale('log')
        ax.set_xlim([2., 10])
        #ax.set_ylim([1e-16, 1e-12])
        #plt.legend(loc = 'best', fontsize = 10)
        ax.grid(True)
        if save:
            print("Saving the figure of the electron Rate fits in ./RateFit_elec.pdf")
            plt.savefig("./RateFit_elec.pdf", bbox_inches='tight')
        else:
            plt.show()
            
    def plotOmegaM(self, save = True):
    
        from matplotlib.lines import Line2D

        # Plot the data
        f, ax = plt.subplots()

        for collision in self.ionData:
            for iMach, Mach in enumerate(collision.Mach):
                alpha_color = 1. - (iMach/(len(collision.Mach)*1.3))**0.5
                ax.plot(collision.Temperature[:],  collision.OmegaM[iMach,:], color=collision.Color, alpha = alpha_color)
            
        #ax.text(0.15, 2e-15, r'$M_{+g} = 15$')#, color=sigma_Xen_ce.Color, alpha = alpha_color+0.2)
        #ax.text(0.25, 1.4e-16, r'$M_{+g} = 0$')#, color=sigma_Xen_ce.Color)

        #custom_lines = [Line2D([0], [0], color=sigma_Xen_iso.Color, lw=2),
        #                Line2D([0], [0], color=sigma_Xen_ce.Color, lw=2)]

        #ax.legend(custom_lines, ['Isotropic', 'Charge-exchange'],loc = 'best', fontsize = 10)
        #plt.title("Interpolation")

        ax.set_ylabel(r'$\Omega^{(1,1,m)}$ $[m^{3}/s]$')
        ax.set_xlabel(r'$T_{+g}$ [eV]')
        ax.set_yscale('log')
        ax.set_xlim([1e-2, 0.4])
        ax.grid(True)
        
        if save:
            print("Saving the figure of the plotOmegaM Rate in ./OmegaM_ion.pdf")
            plt.savefig("./OmegaM_ion.pdf", bbox_inches='tight')
        else:
            plt.show()




