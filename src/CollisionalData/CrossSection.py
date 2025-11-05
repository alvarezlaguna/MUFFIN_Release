import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import random
import colorsys
import numpy as np

import scipy.constants as phy_const
from scipy import integrate
from scipy import interpolate

from scipy.optimize import curve_fit 
import sys
import pickle as pkl
import os.path
np.set_printoptions(threshold=np.inf)


# Function to choose colors of curves
def random_color():
        rand_h = random.randint(0, 100)/100
        rand_sv  = random.randint(90, 100)/100
        return colorsys.hsv_to_rgb(rand_h, rand_sv, rand_sv)

# Class that defines functionalities of cross ssections
class CrossSection:
    
    def __init__(self, name=None, mass=phy_const.m_e, mechanism="e-neutral", potential = 0.0, category=None, density = None, temperature = None,
                 species = None, process = None, param = None, comment = None, updated = None):
        # Attributes
        self.name = name
        self.mass = mass
        self.mechanism = mechanism
        self.potential = potential
        self.category  = category
        self.species   = species
        self.process   = process
        self.param     = param
        self.comment   = comment
        self.updated   = updated
        # Define the energy
        self.Energy = []
        # Define the data
        self.Sigma  = []
        # Define the velocity
        self.Velocity = []
        # Define the temperature for the rate
        if(mechanism == "e-neutral"):
            self.Temperature = np.linspace(1.,20,1000)
            self.Delta = np.linspace(-1, 1, 101)
        if(mechanism == "e-e"):
            self.Delta = np.linspace(-1, 1, 101)
            self.KDeltaZero = np.zeros((101,1000))
            self.KDeltaOne  = np.zeros((101,1000))
            self.NbEnergyPoints = 10000
            self.NbPoints    = 1000
            self.K           = np.zeros(self.NbPoints)
            self.avSigma     = np.zeros(self.NbPoints)
            self.Temperature = np.linspace(1.,20,self.NbPoints)
            self.Energy     = 10**(np.linspace(-3,3,self.NbEnergyPoints))
            self.Velocity   = np.sqrt(2*phy_const.e*self.Energy/self.mass)
            b0Ovg2          = (phy_const.e**2/(4*np.pi*phy_const.epsilon_0*phy_const.m_e/2))
            # Default values that are overriden by the options
            self.T_i        = 4.0  #eV
            self.density    = 1e17 #m^-3
            if temperature is not None:
                self.T_i           = temperature
            if density is not None:
                self.density    = density
            rD2           = phy_const.epsilon_0*self.T_i/(self.density*phy_const.e)
            # Version of the simplified collision with the limiter
            #self.Sigma    = 4*np.pi*(b0Ovg2/self.Velocity**2)**2*np.log(np.sqrt(rD2*self.Velocity**4/b0Ovg2**2))
            #self.Sigma[self.Sigma<0] = 0
            self.Sigma    = 2*np.pi*(b0Ovg2/self.Velocity**2)**2*np.log(1 + (rD2*self.Velocity**4/b0Ovg2**2))
            
        if(mechanism == "ion-neutral"):
            self.Temperature  = np.linspace(1e-2,1,100) 
            self.Mach         = np.linspace(0,15,16)
            self.OmegaM       = np.zeros((16,100))
            self.K_fixedT     = np.zeros(1000)
            self.OmegaE       = np.zeros((16,100))
            self.CoeffsOmegaM = []
            self.CoeffsOmegaE = []
        # Define the rate
        if(mechanism == "e-neutral"):
            self.K_fixedT     = np.zeros(1000)
            self.K      = np.linspace(1.5,20,1000)
            self.KDeltaZero = np.zeros((101,1000))
            self.KDeltaOne  = np.zeros((101,1000))
            self.avSigma     = np.zeros(1000)
        if(mechanism == "ion-neutral"):
            self.K = np.linspace(1e-2,20,100)
        # Define the Omega11
        self.Omega11 = np.linspace(1e-2,20,1000)
        # Define the Omega12
        self.Omega12 = np.linspace(1e-2,20,1000)
        # Define the Omega13
        self.Omega13 = np.linspace(1e-2,20,1000)
        # Define the Omega13
        self.Omega14 = np.linspace(1e-2,20,1000)
        # Define the Omega22
        self.Omega22 = np.linspace(1e-2,20,1000)
        # Define the Omega23
        self.Omega23 = np.linspace(1e-2,20,1000)
        # Define the Omega24
        self.Omega24 = np.linspace(1e-2,20,1000)
        # Define the Omega21
        self.Omega21 = np.linspace(1e-2,20,1000)
        # Define the fit coefficients for the rate
        self.Coeffs = []
        # Define coeffs for sigma
        self.CoeffsSigma = []
        # Define the fit coefficients for the rate
        self.CoeffsOmega11 = []
        # Define the fit coefficients for the rate
        self.CoeffsOmega12 = []
        # Define the fit coefficients for the rate
        self.CoeffsOmega13 = []
        # Define the fit coefficients for the rate
        self.CoeffsOmega14 = []
        # Define the fit coefficients for the rate
        self.CoeffsOmega22 = []
        # Define the fit coefficients for the rate
        self.CoeffsOmega23 = []
        # Define the fit coefficients for the rate
        self.CoeffsOmega24 = []
        # Define the fit coefficients for the rate
        self.CoeffsOmega21 = []
        # Define the color of the curve
        self.Color = random_color()
    
    def readData(self, fileName, units='m2'):
        import csv
        # Automatic detection of delimiter
        delimiter = "\t"
        try:
            with open(fileName, 'r') as myCsvfile:
                header=myCsvfile.readline()
                sniffer = csv.Sniffer()
                dialect = sniffer.sniff(header)
                delimiter = dialect.delimiter
        except Exception:
            pass
        Q = np.loadtxt(fileName, delimiter=delimiter, unpack=True)
        self.Energy   = Q[0]
        if(units=='m2') :
            self.Sigma    = Q[1]
        elif(units=='a2') :
            self.Sigma    = Q[1]*1e-20
        else:
            print("Units not implemented")
            self.Sigma    = Q[1]
        if self.mechanism == "e-neutral":
            self.Velocity = np.sqrt(2*phy_const.e*self.Energy/self.mass)
        if self.mechanism == "ion-neutral":
            self.Velocity = np.sqrt(4*phy_const.e*self.Energy/self.mass)
        self.fileName = fileName
        
    def plotSigma(self, ax):
        ax.plot(self.Energy, self.Sigma, label = self.name, color=self.Color)
        
    def computeRate(self):
        if self.mechanism != "e-e":
            for i, T in enumerate(self.Temperature):
                functional = np.exp(-self.mass*self.Velocity**2/(2*phy_const.e*T))*self.Velocity*self.Sigma*4*np.pi*self.Velocity**2
                self.K[i] = integrate.trapz(functional, x= self.Velocity)*(self.mass/(2*np.pi*phy_const.e*T))**(3/2)
        if self.mechanism == "e-e":
            self.Energy   = 10**(np.linspace(-3,3,self.NbEnergyPoints))
            self.Velocity = np.sqrt(2*phy_const.e*self.Energy/self.mass)
            b0Ovg2        = (phy_const.e**2/(4*np.pi*phy_const.epsilon_0*phy_const.m_e/2))
            for i, T in enumerate(self.Temperature):
                T_e           = T
                n_e           = self.density
                rD2           = phy_const.epsilon_0*T_e/(n_e*phy_const.e)
                self.Sigma    = 2*np.pi*(b0Ovg2/self.Velocity**2)**2*np.log(1 + (rD2*self.Velocity**4/b0Ovg2**2))
                Beta_e_Half   = phy_const.m_e/(2*phy_const.e*T_e)/2
                functional = np.exp(-Beta_e_Half*self.Velocity**2)*self.Velocity*self.Sigma*4*np.pi*self.Velocity**2
                self.K[i] = integrate.trapz(functional, x= self.Velocity)*(1/np.pi*Beta_e_Half)**(3/2)
                
    def computeRateDeltaZero(self):
        if self.mechanism != "e-e":
            def computeKInelasticApprox(T_e, Delta_e, order):

                m_e     = phy_const.m_e
                e       = phy_const.e
                beta_e  = m_e/(2*e*T_e)

                # Definitions of the Electrons VDF
                a_0 = (1 + 15./8.*Delta_e)
                a_2 = -5./2.*beta_e*Delta_e
                a_4 = beta_e*beta_e/2.*Delta_e

                vel = self.Velocity

                f_eOvNe = (beta_e/np.pi)**(3./2.)*np.exp(-beta_e*vel**2)*(a_0 + a_2*vel**2 + a_4*vel**4)

                # Limit to the last positive value
                lastPositive = 0
                for idx, ifunc in enumerate(f_eOvNe):
                    check = True
                    if ifunc > 0 and check:
                        f_eOvNe[idx] = ifunc
                        lastPositive = idx

                    else:
                        f_eOvNe[idx] = 0
                        check = False

                functional = vel[0:lastPositive]**(2*order + 3)*f_eOvNe[0:lastPositive]*self.Sigma[0:lastPositive]

                freq = 4*np.pi*integrate.trapz(functional, x= vel[0:lastPositive])*(2*beta_e)**order

                return freq
            
            fileName = "./"+os.path.splitext(self.fileName)[0]+"_KDelta_0.pkl"
            if os.path.exists(fileName):
                print("reading the rates for ", self.name," in ",fileName)
                with open(fileName,'rb') as f:
                    self.KDeltaZero = pkl.load(f)
            else:
                for iDelta, Del in enumerate(self.Delta):
                    print("Computing Delta = ", Del)
                    for iTemp, T in enumerate(self.Temperature):
                        self.KDeltaZero[iDelta, iTemp] = computeKInelasticApprox(T, Del, 0)
                print("saving the rates for ", self.name," in ",fileName)
                with open(fileName,'wb') as f:
                    pkl.dump(self.KDeltaZero, f)
                
    def computeRateDeltaOne(self):
        if self.mechanism != "e-e":
            def computeKInelasticApprox(T_e, Delta_e, order):

                m_e     = phy_const.m_e
                e       = phy_const.e
                beta_e  = m_e/(2*e*T_e)

                # Definitions of the Electrons VDF
                a_0 = (1 + 15./8.*Delta_e)
                a_2 = -5./2.*beta_e*Delta_e
                a_4 = beta_e*beta_e/2.*Delta_e

                vel = self.Velocity

                f_eOvNe = (beta_e/np.pi)**(3./2.)*np.exp(-beta_e*vel**2)*(a_0 + a_2*vel**2 + a_4*vel**4)

                # Limit to the last positive value
                lastPositive = 0
                for idx, ifunc in enumerate(f_eOvNe):
                    check = True
                    if ifunc > 0 and check:
                        f_eOvNe[idx] = ifunc
                        lastPositive = idx

                    else:
                        f_eOvNe[idx] = 0
                        check = False

                functional = vel[0:lastPositive]**(2*order + 3)*f_eOvNe[0:lastPositive]*self.Sigma[0:lastPositive]

                freq = 4*np.pi*integrate.trapz(functional, x= vel[0:lastPositive])*(2*beta_e)**order

                return freq
            
            fileName = "./"+os.path.splitext(self.fileName)[0]+"_KDelta_1.pkl"
            if os.path.exists(fileName):
                print("reading the rates for ", self.name," in ",fileName)
                with open(fileName,'rb') as f:
                    self.KDeltaOne = pkl.load(f)
            else:
                for iDelta, Del in enumerate(self.Delta):
                    print("Computing Delta = ", Del)
                    for iTemp, T in enumerate(self.Temperature):
                        self.KDeltaOne[iDelta, iTemp] = computeKInelasticApprox(T, Del, 1)
                print("saving the rates for ", self.name," in ",fileName)
                with open(fileName,'wb') as f:
                    pkl.dump(self.KDeltaOne, f)
                
    def computeKInelasticApprox(self, T_e, Delta_e, order):

        m_e     = phy_const.m_e
        e       = phy_const.e
        beta_e  = m_e/(2*e*T_e)

        # Definitions of the Electrons VDF
        a_0 = (1 + 15./8.*Delta_e)
        a_2 = -5./2.*beta_e*Delta_e
        a_4 = beta_e*beta_e/2.*Delta_e

        vel = self.Velocity

        f_eOvNe = (beta_e/np.pi)**(3./2.)*np.exp(-beta_e*vel**2)*(a_0 + a_2*vel**2 + a_4*vel**4)

        # Limit to the last positive value
        lastPositive = 0
        for idx, ifunc in enumerate(f_eOvNe):
            check = True
            if ifunc > 0 and check:
                f_eOvNe[idx] = ifunc
                lastPositive = idx

            else:
                f_eOvNe[idx] = 0
                check = False

        functional = vel[0:lastPositive]**(2*order + 3)*f_eOvNe[0:lastPositive]*self.Sigma[0:lastPositive]

        freq = 4*np.pi*integrate.trapz(functional, x= vel[0:lastPositive])*(2*beta_e)**order

        return freq


                
    def computeAverageSigma(self):
        if self.mechanism != "e-e":
            for i, T in enumerate(self.Temperature):
                functional = np.exp(-self.mass*self.Velocity**2/(2*phy_const.e*T))*self.Sigma*4*np.pi*self.Velocity**2
                self.avSigma[i] = integrate.trapz(functional, x= self.Velocity)*(self.mass/(2*np.pi*phy_const.e*T))**(3/2)
        if self.mechanism == "e-e":
            for i, T in enumerate(self.Temperature):
                self.Energy   = 10**(np.linspace(-3,3,self.NbEnergyPoints))
                self.Velocity = np.sqrt(2*phy_const.e*self.Energy/self.mass)
                b0Ovg2        = (phy_const.e**2/(4*np.pi*phy_const.epsilon_0*phy_const.m_e/2))
                T_e           = T
                n_e           = self.density
                rD2           = phy_const.epsilon_0*T_e/(n_e*phy_const.e)
                self.Sigma    = 2*np.pi*(b0Ovg2/self.Velocity**2)**2*np.log(1 + (rD2*self.Velocity**4/b0Ovg2**2))
                functional = np.exp(-self.mass*self.Velocity**2/(2*phy_const.e*T))*self.Sigma*4*np.pi*self.Velocity**2
                self.avSigma[i] = integrate.trapz(functional, x= self.Velocity)*(self.mass/(2*np.pi*phy_const.e*T))**(3/2)
    
    def computeRateFunction(self):
        if self.mechanism != "e-e":
            functional = np.exp(-self.mass*self.Velocity**2/(2*phy_const.e*temperature))*self.Velocity*self.Sigma*4*np.pi*self.Velocity**2
            rate = integrate.trapz(functional, x= self.Velocity)*(self.mass/(2*np.pi*phy_const.e*temperature))**(3/2)
        if self.mechanism == "e-e":
            self.Energy   = 10**(np.linspace(-3,3,self.NbEnergyPoints))
            self.Velocity = np.sqrt(2*phy_const.e*self.Energy/self.mass)
            b0Ovg2        = (phy_const.e**2/(4*np.pi*phy_const.epsilon_0*phy_const.m_e/2))
            T_e           = temperature
            n_e           = self.density
            rD2           = phy_const.epsilon_0*T_e/(n_e*phy_const.e)
            self.Sigma    = 2*np.pi*(b0Ovg2/self.Velocity**2)**2*np.log(1 + (rD2*self.Velocity**4/b0Ovg2**2))
            functional    = np.exp(-self.mass*self.Velocity**2/(2*phy_const.e*T))*self.Velocity*self.Sigma*4*np.pi*self.Velocity**2
            rate = integrate.trapz(functional, x= self.Velocity)*(self.mass/(2*np.pi*phy_const.e*T))**(3/2)
        return rate
        
        
    def plotRate(self, ax):
        ax.plot(self.Temperature, self.K, label = self.name, color=self.Color)
        
    def fitRate(self):
        if(self.mechanism == "e-neutral" or self.mechanism == "e-e"):
            Approx = curve_fit(lambda t,a,b,c,d,f: a*t**f*np.exp(-b*t + c*t**2 - d*t**3.),  1./self.Temperature,  self.K,  p0=(self.K[-1], 0.,0.,0.,0.))
            self.Coeffs = curve_fit(lambda t,a,b,c,d,e,f: a*t**f*np.exp(-b*t+ c*t**2 - d*t**3. + e*t**4),  1./self.Temperature,  self.K,  p0=(Approx[0][0], Approx[0][1], Approx[0][2], Approx[0][3], 0., Approx[0][-1]), method='dogbox')[0]
        if(self.mechanism == "ion-neutral"):
            self.Coeffs = curve_fit(lambda t,a,b,c,d: a*t**d*np.exp(-b*t+ c*t**2),  1./self.Temperature,  self.K,  p0=(self.K[-1], 0.01, 0.01, 0.01))[0]

    def plotFitRate(self, ax):
        def fit_function(T):
            if(self.mechanism == "e-neutral" or self.mechanism == "e-e"):
                a = self.Coeffs[0]
                b = self.Coeffs[1]
                c = self.Coeffs[2]
                d = self.Coeffs[3]
                e = self.Coeffs[4]
                f = self.Coeffs[5]

                t = (1/T)

                return a*t**f*np.exp(-b*t+ c*t**2 - d*t**3. + e*t**4)
            if(self.mechanism == "ion-neutral"):
                a = self.Coeffs[0]
                b = self.Coeffs[1]
                c = self.Coeffs[2]
                d = self.Coeffs[3]

                t = (1/T)

                return a*t**d*np.exp(-b*t+ c*t**2)
        
        fit = fit_function(self.Temperature)
        ax.plot(self.Temperature, fit, label = self.name, linestyle='None', marker='x', markevery=15, color=self.Color)

        
    def computeOmega11(self):
        for i, T in enumerate(self.Temperature):
            psi             = np.sqrt(self.mass*self.Velocity**2/(2*phy_const.e*T))
            Sigma           = self.Sigma/(2*np.pi)
            functional      = psi**5*np.exp(-psi**2)*Sigma
            self.Omega11[i] = integrate.trapz(functional, x= psi)*np.sqrt(2*np.pi*phy_const.e*T/self.mass)

    def computeOmega12(self):
        for i, T in enumerate(self.Temperature):
            psi             = np.sqrt(self.mass*self.Velocity**2/(2*phy_const.e*T))
            Sigma           = self.Sigma/(2*np.pi)
            functional      = psi**7*np.exp(-psi**2)*Sigma
            self.Omega12[i] = integrate.trapz(functional, x= psi)*np.sqrt(2*np.pi*phy_const.e*T/self.mass)

    def computeOmega13(self):
        for i, T in enumerate(self.Temperature):
            psi             = np.sqrt(self.mass*self.Velocity**2/(2*phy_const.e*T))
            Sigma           = self.Sigma/(2*np.pi)
            functional      = psi**9*np.exp(-psi**2)*Sigma
            self.Omega13[i] = integrate.trapz(functional, x= psi)*np.sqrt(2*np.pi*phy_const.e*T/self.mass)
            
    def computeOmega14(self):
        for i, T in enumerate(self.Temperature):
            psi             = np.sqrt(self.mass*self.Velocity**2/(2*phy_const.e*T))
            Sigma           = self.Sigma/(2*np.pi)
            functional      = psi**11*np.exp(-psi**2)*Sigma
            self.Omega14[i] = integrate.trapz(functional, x= psi)*np.sqrt(2*np.pi*phy_const.e*T/self.mass)

    def computeOmega22(self):
        if self.mechanism == "e-e":
            for i, T in enumerate(self.Temperature):
                    Beta_e_Half   = self.mass/(2*phy_const.e*T)/2.
                    psi           = np.sqrt(self.Velocity**2*Beta_e_Half)
                    b0Ovg2        = (phy_const.e**2/(4*np.pi*phy_const.epsilon_0*phy_const.m_e/2.))
                    T_e           = T
                    n_e           = self.density
                    rD2           = phy_const.epsilon_0*T_e/(n_e*phy_const.e)
                    
                    # Simplified version
                    #Sigma2        = 8*np.pi*(b0Ovg2/self.Velocity**2)**2*np.log(np.sqrt(rD2*self.Velocity**4/b0Ovg2**2))
                    # Full Version
                    #Sigma2        = 4*np.pi*(b0Ovg2/self.Velocity**2)**2*(np.log(1 + (rD2*self.Velocity**4/b0Ovg2**2)) - (rD2*self.Velocity**4/b0Ovg2**2)/(1 + (rD2*self.Velocity**4/b0Ovg2**2)))
                    #Sigma2        = 4*np.pi*(b0Ovg2/self.Velocity**2)**2*(np.log(1 + (rD2*self.Velocity**4/b0Ovg2**2)) - (rD2*self.Velocity**4/b0Ovg2**2)/(1 + (rD2*self.Velocity**4/b0Ovg2**2)))
                    Lambda          = 4*np.pi*phy_const.epsilon_0/phy_const.e**2*phy_const.m_e/2.*np.sqrt(rD2)*3/(Beta_e_Half/2.)
                    Sigma2          = 8*np.pi*(b0Ovg2/self.Velocity**2)**2*(np.log(Lambda))# + np.log(2.) - 1. - 0.5772)
                    mu_ee           = self.mass/2.
                    functional      = psi**7*np.exp(-psi**2)*Sigma2
                    self.Omega22[i] = 1./2.*integrate.trapz(functional, x= psi)*np.sqrt(1/np.pi*1/Beta_e_Half)
        else:
            for i, T in enumerate(self.Temperature):
                psi             = np.sqrt(self.mass*self.Velocity**2/(2*phy_const.e*T))
                Sigma2          = 2./3.*self.Sigma/(2*np.pi)
                functional      = psi**7*np.exp(-psi**2)*Sigma2
                self.Omega22[i] = integrate.trapz(functional, x= psi)*np.sqrt(2*np.pi*phy_const.e*T/self.mass)
    
    def computeOmega23(self):
        if self.mechanism == "e-e":
            for i, T in enumerate(self.Temperature):
                    Beta_e_Half   = self.mass/(2*phy_const.e*T)/2
                    psi           = np.sqrt(self.Velocity**2*Beta_e_Half)
                    b0Ovg2        = (phy_const.e**2/(4*np.pi*phy_const.epsilon_0*phy_const.m_e/2))
                    T_e           = T
                    n_e           = self.density
                    rD2           = phy_const.epsilon_0*T_e/(n_e*phy_const.e)
                    
                    # Simplified version
                    #Sigma2        = 8*np.pi*(b0Ovg2/self.Velocity**2)**2*np.log(np.sqrt(rD2*self.Velocity**4/b0Ovg2**2))
                    # Full Version
                    Sigma2        = 4*np.pi*(b0Ovg2/self.Velocity**2)**2*(np.log(1 + (rD2*self.Velocity**4/b0Ovg2**2)) - (rD2*self.Velocity**4/b0Ovg2**2)/(1 + (rD2*self.Velocity**4/b0Ovg2**2)))
                    #Sigma2        = 4*np.pi*(b0Ovg2/self.Velocity**2)**2*(np.log(1 + (rD2*self.Velocity**4/b0Ovg2**2)) - (rD2*self.Velocity**4/b0Ovg2**2)/(1 + (rD2*self.Velocity**4/b0Ovg2**2)))
                    mu_ee           = self.mass/2.
                    functional      = psi**(2*3 + 3)*np.exp(-psi**2)*Sigma2
                    self.Omega23[i] = 1./2.*integrate.trapz(functional, x= psi)*np.sqrt(1/np.pi*1/Beta_e_Half)
        else:
            for i, T in enumerate(self.Temperature):
                psi             = np.sqrt(self.mass*self.Velocity**2/(2*phy_const.e*T))
                Sigma2          = 2./3.*self.Sigma/(2*np.pi)
                functional      = psi**(2*3 + 3)*np.exp(-psi**2)*Sigma2
                self.Omega23[i] = integrate.trapz(functional, x= psi)*np.sqrt(2*np.pi*phy_const.e*T/self.mass)
    
    def computeOmega24(self):
        if self.mechanism == "e-e":
            for i, T in enumerate(self.Temperature):
                    Beta_e_Half   = self.mass/(2*phy_const.e*T)/2
                    psi           = np.sqrt(self.Velocity**2*Beta_e_Half)
                    b0Ovg2        = (phy_const.e**2/(4*np.pi*phy_const.epsilon_0*phy_const.m_e/2))
                    T_e           = T
                    n_e           = self.density
                    rD2           = phy_const.epsilon_0*T_e/(n_e*phy_const.e)
                    
                    # Simplified version
                    #Sigma2        = 8*np.pi*(b0Ovg2/self.Velocity**2)**2*np.log(np.sqrt(rD2*self.Velocity**4/b0Ovg2**2))
                    # Full Version
                    Sigma2        = 4*np.pi*(b0Ovg2/self.Velocity**2)**2*(np.log(1 + (rD2*self.Velocity**4/b0Ovg2**2)) - (rD2*self.Velocity**4/b0Ovg2**2)/(1 + (rD2*self.Velocity**4/b0Ovg2**2)))
                    #Sigma2        = 4*np.pi*(b0Ovg2/self.Velocity**2)**2*(np.log(1 + (rD2*self.Velocity**4/b0Ovg2**2)) - (rD2*self.Velocity**4/b0Ovg2**2)/(1 + (rD2*self.Velocity**4/b0Ovg2**2)))
                    mu_ee           = self.mass/2.
                    functional      = psi**(2*4 + 3)*np.exp(-psi**2)*Sigma2
                    self.Omega24[i] = 1./2.*integrate.trapz(functional, x= psi)*np.sqrt(1/np.pi*1/Beta_e_Half)
        else:
            for i, T in enumerate(self.Temperature):
                psi             = np.sqrt(self.mass*self.Velocity**2/(2*phy_const.e*T))
                Sigma2          = 2./3.*self.Sigma/(2*np.pi)
                functional      = psi**(2*4 + 3)*np.exp(-psi**2)*Sigma2
                self.Omega24[i] = integrate.trapz(functional, x= psi)*np.sqrt(2*np.pi*phy_const.e*T/self.mass)
        

    def computeOmega21(self):
        for i, T in enumerate(self.Temperature):
            psi             = np.sqrt(self.mass*self.Velocity**2/(2*phy_const.e*T))
            Sigma2          = 2./3.*self.Sigma/(2*np.pi)
            functional      = psi**5*np.exp(-psi**2)*Sigma2
            self.Omega21[i] = integrate.trapz(functional, x= psi)*np.sqrt(2*np.pi*phy_const.e*T/self.mass)

    def fitOmega11(self):
        self.CoeffsOmega11 = curve_fit(lambda t,a,b,c,d,e,f: a*t**f*np.exp(-b*t+ c*t**2 - d*t**3 + e*t**4),  1./self.Temperature,  self.Omega11,  p0=(self.Omega11[-1], 0.01, 0.01, 0.01, 0.01, 0.001))
    
    def fitOmega12(self):
        self.CoeffsOmega12 = curve_fit(lambda t,a,b,c,d,e,f: a*t**f*np.exp(-b*t+ c*t**2 - d*t**3 + e*t**4),  1./self.Temperature,  self.Omega12,  p0=(self.Omega12[-1], 0.01, 0.01, 0.01, 0.01, 0.001))
    
    def fitOmega13(self):
        self.CoeffsOmega13 = curve_fit(lambda t,a,b,c,d,e,f: a*t**f*np.exp(-b*t+ c*t**2 - d*t**3 + e*t**4),  1./self.Temperature,  self.Omega13,  p0=(self.Omega13[-1], 0.01, 0.01, 0.01, 0.01, 0.001))
    
    def fitOmega14(self):
        self.CoeffsOmega14 = curve_fit(lambda t,a,b,c,d,e,f: a*t**f*np.exp(-b*t+ c*t**2 - d*t**3 + e*t**4),  1./self.Temperature,  self.Omega14,  p0=(self.Omega14[-1], 0.01, 0.01, 0.01, 0.01, 0.001))
    
    def fitOmega22(self):
        self.CoeffsOmega22 = curve_fit(lambda t,a,b,c,d,e,f: a*t**f*np.exp(-b*t+ c*t**2 - d*t**3 + e*t**4),  1./self.Temperature,  self.Omega22,  p0=(self.Omega22[-1], 0.01, 0.01, 0.01, 0.01, 0.001))
    def fitOmega23(self):
        self.CoeffsOmega23 = curve_fit(lambda t,a,b,c,d,e,f: a*t**f*np.exp(-b*t+ c*t**2 - d*t**3 + e*t**4),  1./self.Temperature,  self.Omega23,  p0=(self.Omega23[-1], 0.01, 0.01, 0.01, 0.01, 0.001))
    def fitOmega24(self):
        self.CoeffsOmega24 = curve_fit(lambda t,a,b,c,d,e,f: a*t**f*np.exp(-b*t+ c*t**2 - d*t**3 + e*t**4),  1./self.Temperature,  self.Omega24,  p0=(self.Omega24[-1], 0.01, 0.01, 0.01, 0.01, 0.001))
    
    def fitOmega21(self):
        self.CoeffsOmega21 = curve_fit(lambda t,a,b,c,d,e,f: a*t**f*np.exp(-b*t+ c*t**2 - d*t**3 + e*t**4),  1./self.Temperature,  self.Omega21,  p0=(self.Omega21[-1], 0.01, 0.01, 0.01, 0.01, 0.001))
    
    
    def plotFitOmega11(self, ax):
        def fit_function(T):
            a = self.CoeffsOmega11[0][0]
            b = self.CoeffsOmega11[0][1]
            c = self.CoeffsOmega11[0][2]
            d = self.CoeffsOmega11[0][3]
            e = self.CoeffsOmega11[0][4]
            f = self.CoeffsOmega11[0][5]
            
            t = (1/T)
            
            return a*t**f*np.exp(-b*t+ c*t**2 - d*t**3 + e*t**4)
        
        fit = fit_function(self.Temperature)
        ax.plot(self.Temperature, fit, label = self.name, linestyle='None', marker='x', markevery=15, color=self.Color)
    
    def plotFitOmega12(self, ax):
        def fit_function(T):
            a = self.CoeffsOmega12[0][0]
            b = self.CoeffsOmega12[0][1]
            c = self.CoeffsOmega12[0][2]
            d = self.CoeffsOmega12[0][3]
            e = self.CoeffsOmega12[0][4]
            f = self.CoeffsOmega12[0][5]
            
            t = (1/T)
            
            return a*t**f*np.exp(-b*t+ c*t**2 - d*t**3 + e*t**4)
        
        fit = fit_function(self.Temperature)
        ax.plot(self.Temperature, fit, label = self.name, linestyle='None', marker='x', markevery=15, color=self.Color)
    
    
    def plotFitOmega13(self, ax):
        def fit_function(T):
            a = self.CoeffsOmega13[0][0]
            b = self.CoeffsOmega13[0][1]
            c = self.CoeffsOmega13[0][2]
            d = self.CoeffsOmega13[0][3]
            e = self.CoeffsOmega13[0][4]
            f = self.CoeffsOmega13[0][5]
            
            t = (1/T)
            
            return a*t**f*np.exp(-b*t+ c*t**2 - d*t**3 + e*t**4)
        
        fit = fit_function(self.Temperature)
        ax.plot(self.Temperature, fit, label = self.name, linestyle='None', marker='x', markevery=15, color=self.Color)

    def plotFitOmega14(self, ax):
       def fit_function(T):
           a = self.CoeffsOmega14[0][0]
           b = self.CoeffsOmega14[0][1]
           c = self.CoeffsOmega14[0][2]
           d = self.CoeffsOmega14[0][3]
           e = self.CoeffsOmega14[0][4]
           f = self.CoeffsOmega14[0][5]
           
           t = (1/T)
           
           return a*t**f*np.exp(-b*t+ c*t**2 - d*t**3 + e*t**4)
       
       fit = fit_function(self.Temperature)
       ax.plot(self.Temperature, fit, label = self.name, linestyle='None', marker='x', markevery=15, color=self.Color)

    
    def plotFitOmega22(self, ax):
        def fit_function(T):
            a = self.CoeffsOmega22[0][0]
            b = self.CoeffsOmega22[0][1]
            c = self.CoeffsOmega22[0][2]
            d = self.CoeffsOmega22[0][3]
            e = self.CoeffsOmega22[0][4]
            f = self.CoeffsOmega22[0][5]
            
            t = (1/T)
            
            return a*t**f*np.exp(-b*t+ c*t**2 - d*t**3 + e*t**4)
        
        fit = fit_function(self.Temperature)
        ax.plot(self.Temperature, fit, label = self.name, linestyle='None', marker='x', markevery=15, color=self.Color)
    
    def plotFitOmega23(self, ax):
        def fit_function(T):
            a = self.CoeffsOmega23[0][0]
            b = self.CoeffsOmega23[0][1]
            c = self.CoeffsOmega23[0][2]
            d = self.CoeffsOmega23[0][3]
            e = self.CoeffsOmega23[0][4]
            f = self.CoeffsOmega23[0][5]
            
            t = (1/T)
            
            return a*t**f*np.exp(-b*t+ c*t**2 - d*t**3 + e*t**4)
        
        fit = fit_function(self.Temperature)
        ax.plot(self.Temperature, fit, label = self.name, linestyle='None', marker='x', markevery=15, color=self.Color)
    
    def plotFitOmega24(self, ax):
        def fit_function(T):
            a = self.CoeffsOmega24[0][0]
            b = self.CoeffsOmega24[0][1]
            c = self.CoeffsOmega24[0][2]
            d = self.CoeffsOmega24[0][3]
            e = self.CoeffsOmega24[0][4]
            f = self.CoeffsOmega24[0][5]
            
            t = (1/T)
            
            return a*t**f*np.exp(-b*t+ c*t**2 - d*t**3 + e*t**4)
        
        fit = fit_function(self.Temperature)
        ax.plot(self.Temperature, fit, label = self.name, linestyle='None', marker='x', markevery=15, color=self.Color)
 

    
    def plotOmega11(self, ax):
        ax.plot(self.Temperature, self.Omega11, label = self.name, color=self.Color)
    
    def plotOmega12(self, ax):
        ax.plot(self.Temperature, self.Omega12, label = self.name, color=self.Color)
    
    def plotOmega13(self, ax):
        ax.plot(self.Temperature, self.Omega13, label = self.name, color=self.Color)
    
    def plotOmega14(self, ax):
        ax.plot(self.Temperature, self.Omega14, label = self.name, color=self.Color)
    
    def plotOmega22(self, ax):
        ax.plot(self.Temperature, self.Omega22, label = self.name, color=self.Color)

    def plotOmega23(self, ax):
        ax.plot(self.Temperature, self.Omega23, label = self.name, color=self.Color)

    def plotOmega24(self, ax):
        ax.plot(self.Temperature, self.Omega24, label = self.name, color=self.Color)

    def plotNu22(self, ax):
        ax.plot(self.Temperature, 16./5.*np.sqrt((np.pi*phy_const.e*self.Temperature)/self.mass)*self.Omega22, label = self.name, color=self.Color, linestyle='--')


    def computeOmegaM(self):
        fileName = "./"+os.path.splitext(self.fileName)[0]+"_OmegaM.pkl"
        if os.path.exists(fileName):
            print("reading the rates for ", self.name," in ",fileName)
            with open(fileName,'rb') as f:
                self.OmegaM = pkl.load(f)
        else:
            for j, Mach in enumerate(self.Mach):
                print("Computing Mach = ", Mach)
                for i, T in enumerate(self.Temperature):
                    
                    m_ig   = 0.5*self.mass
                    psi    = np.sqrt(m_ig*self.Velocity**2/(2*phy_const.e*T))
                    
                    logPsilarge    = np.linspace(np.log10(1e-12),np.log10(psi[-1]), 10000)
                    logSigma       = np.log10(self.Sigma)
                    
                    logPsi         = np.log10(psi)
                    logPsi[0]      = -12
                    
                    f = interpolate.interp1d(logPsi, logSigma)
                    
                    psiLarge   = 10**(logPsilarge)
                    sigmaLarge = 10**(f(logPsilarge))

                    y = 2*Mach*psiLarge
                    functional = np.zeros(np.shape(psiLarge))
                    
                    for k, i_y in enumerate(y):
                        if Mach == 0 :
                            F_m           = 1 + i_y**2/10.
                            functional[k] = F_m*psiLarge[k]**5*np.exp(-psiLarge[k]**2)*sigmaLarge[k]
                            
                        else:
                            psi_k = psiLarge[k]
                            
                            exponent_1 = -(psi_k - Mach)**2
                            exponent_2 = -(psi_k + Mach)**2

                            firstPart  = (np.exp(exponent_1) + np.exp(exponent_2))
                            secondPart =  (np.exp(exponent_1) - np.exp(exponent_2))

                            functional[k] = psi_k**2*sigmaLarge[k]*3*(firstPart*i_y - secondPart)/(2*(2*Mach)**3)

                    self.OmegaM[j,i] = integrate.trapz(functional, x= psiLarge)
                    self.OmegaM[j,i] = self.OmegaM[j,i]*(phy_const.e*T/(2*np.pi*m_ig))**(1./2.)
            print("saving the rates for ", self.name," in ",fileName)
            with open(fileName,'wb') as f:
                pkl.dump(self.OmegaM, f)
        
    def computeOmegaMBad(self):
        #functional = np.zeros(np.shape(self.Velocity))
        for j, Mach in enumerate(self.Mach):
            print("Computing Mach = ", Mach)
            for i, T in enumerate(self.Temperature):
                
                m_ig   = 0.5*self.mass
                psi    = np.sqrt(m_ig*self.Velocity**2/(2*phy_const.e*T))
                
                logPsilarge    = np.linspace(np.log10(1e-12),np.log10(psi[-1]), 100)
                logSigma       = np.log10(self.Sigma)
                
                logPsi         = np.log10(psi)
                logPsi[0]      = -12
                
                f = interpolate.interp1d(logPsi, logSigma)
                
                psiLarge   = 10**(logPsilarge)
                sigmaLarge = 10**(f(logPsilarge))

                y = 2*Mach*psiLarge
                functional = np.zeros(np.shape(psiLarge))
                
                for k, i_y in enumerate(y):
                    if Mach == 0 : 
                        F_m           = 1 + i_y**2/10.
                        functional[k] = F_m*psiLarge[k]**5*np.exp(-psiLarge[k]**2)*sigmaLarge[k]
                        
                    else:
                        psi_k = psiLarge[k]
                        
                        exponent_1 = -(psi_k - Mach)**2
                        exponent_2 = -(psi_k + Mach)**2

                        firstPart  = (np.exp(exponent_1) + np.exp(exponent_2))
                        secondPart =  (np.exp(exponent_1) - np.exp(exponent_2))

                        functional[k] = psi_k**2*sigmaLarge[k]*3*(firstPart*i_y - secondPart)/(2*(2*Mach)**3)



                self.OmegaM[j,i] = integrate.trapz(functional, x= psiLarge)
                self.OmegaM[j,i] = self.OmegaM[j,i]*(phy_const.e*T/(2*np.pi*m_ig))**(1./2.) 
    
    def fitOmegaM(self):
        fileName = "./"+os.path.splitext(self.fileName)[0]+"_fitOmegaM.pkl"
        if os.path.exists(fileName):
            print("reading the coefficients for ", self.name,"in ",fileName)
            with open(fileName,'rb') as f:
                self.Mach, self.CoeffsOmegaM = pkl.load(f)
        else:
            for i, Mach in enumerate(self.Mach):
                CoeffsOmegaM_Mach = curve_fit(lambda t,a,b,c,d: a*t**d*np.exp(-b*t+ c*t**2),  1./self.Temperature,  self.OmegaM[i,:],  p0=(self.OmegaM[i,-1], 0.01, 0.01, 0.01))
                self.CoeffsOmegaM.append(CoeffsOmegaM_Mach[0])
            print("saving the coefficients for ", self.name,"in ",fileName)
            with open(fileName,'wb') as f:
                pkl.dump([self.Mach, self.CoeffsOmegaM], f)
                #print(self.CoeffsOmegaM)
        
    
    def computeOmegaE(self):
        fileName = "./"+os.path.splitext(self.fileName)[0]+"_OmegaE.pkl"
        if os.path.exists(fileName):
            print("reading the rates for ", self.name," in ",fileName)
            with open(fileName,'rb') as f:
                self.OmegaE = pkl.load(f)
        else:
            #functional = np.zeros(np.shape(self.Velocity))
            for j, Mach in enumerate(self.Mach):
                print("Computing Mach = ", Mach)
                for i, T in enumerate(self.Temperature):
                    
                    m_ig   = 0.5*self.mass
                    psi    = np.sqrt(m_ig*self.Velocity**2/(2*phy_const.e*T))
                    
                    logPsilarge    = np.linspace(np.log10(1e-12),np.log10(psi[-1]), 10000)
                    logSigma       = np.log10(self.Sigma)
                    
                    logPsi         = np.log10(psi)
                    logPsi[0]      = -12
                    
                    f = interpolate.interp1d(logPsi, logSigma)
                    
                    psiLarge   = 10**(logPsilarge)
                    sigmaLarge = 10**(f(logPsilarge))

                    y = 2*Mach*psiLarge
                    functional = np.zeros(np.shape(psiLarge))
                    
                    for k, i_y in enumerate(y):
                        if Mach == 0 :
                            F_m           = 1 + i_y**2/6.
                            functional[k] = F_m*psiLarge[k]**5*np.exp(-psiLarge[k]**2)*sigmaLarge[k]
                            
                        else:
                            psi_k = psiLarge[k]
                            
                            exponent_1 = -(psi_k - Mach)**2
                            exponent_2 = -(psi_k + Mach)**2

                            firstPart  = (np.exp(exponent_1) - np.exp(exponent_2))

                            functional[k] = psi_k**4*sigmaLarge[k]*firstPart/(2*(2*Mach))

                    self.OmegaE[j,i] = integrate.trapz(functional, x= psiLarge)
                    self.OmegaE[j,i] = self.OmegaE[j,i]*(phy_const.e*T/(2*np.pi*m_ig))**(1./2.)
            
            print("saving the rates for ", self.name," in ",fileName)
            with open(fileName,'wb') as f:
                pkl.dump(self.OmegaE, f)
                
    def fitOmegaE(self):
        fileName = "./"+os.path.splitext(self.fileName)[0]+"_fitOmegaE.pkl"
        if os.path.exists(fileName):
            print("reading the coefficients for ", self.name,"in ",fileName)
            with open(fileName,'rb') as f:
                self.Mach, self.CoeffsOmegaE = pkl.load(f)
        else:
            for i, Mach in enumerate(self.Mach):
                CoeffsOmegaE_Mach = curve_fit(lambda t,a,b,c,d: a*t**d*np.exp(-b*t+ c*t**2),  1./self.Temperature,  self.OmegaE[i,:],  p0=(self.OmegaE[i,-1], 0.01, 0.01, 0.01))
                self.CoeffsOmegaE.append(CoeffsOmegaE_Mach[0])
            print("saving the coefficients for ", self.name,"in ",fileName)
            with open(fileName,'wb') as f:
                pkl.dump([self.Mach, self.CoeffsOmegaE], f)
                #print(self.CoeffsOmegaE)




    def computeK_fixed_T(self, T_s, sigma_cst=False, save=False):
        T_g = 0.025 #eV
        m_g = self.mass + phy_const.m_e
        massRatio = self.mass/m_g
        m_sg   = self.mass/(1 + massRatio)
#        T_g = 0.025
        T_sg   = (T_s+T_g*massRatio)/(1+massRatio) #!!!
        psi    = np.sqrt(m_sg/(phy_const.e*T_sg))*self.Velocity
        logMach = np.linspace(np.log10(psi[1]), np.log10(psi[-1]/5) , 999)
        self.Mach = np.array([0] + list(10**logMach))                       #TODO: change Mach into velocity (?)
        if sigma_cst == True:
            fileName = "./" + os.path.splitext(self.fileName)[0] + "_K_T=%s_Tg=%s_sigma_cst.txt"%(T_s,T_g)
        else:
            #fileName = "./"+os.path.splitext(self.fileName)[0]+"_K_T=%s.pkl"%Ti
            fileName = "./"+os.path.splitext(self.fileName)[0]+"_K_T=%s_Tg=%s.txt"%(T_s, T_g)
            #fileName = "./"+os.path.splitext(self.fileName)[0]+"_K_Ti=%s_Tg=%s.pkl"%(Ti,Tg)
        if os.path.exists(fileName):
            print("reading the rates for ", self.name," in ",fileName)
#            with open(fileName,'rb') as file:
#                self.K_fixedT = pkl.load(file)
            self.K_fixedT = np.loadtxt(fileName)
        else:
            print('computing the rates for ', self.name)
#            logPsilarge    = np.linspace(np.log10(1e-12),np.log10(psi[-1]), 10000)
            logPsilarge    = np.linspace(np.log10(psi[1]), np.log10(psi[-1]), 1000)
            logSigma       = np.log10(self.Sigma)
            logPsi         = np.log10(psi)
            logPsi[0]      = -12
            f = interpolate.interp1d(logPsi, logSigma)#, kind='quadratic')
            for i, Mach in enumerate(self.Mach):
#                print("Computing Mach = ", Mach)
                if Mach > 30:
                    self.K_fixedT[i] = m_sg/self.mass*10**f(np.log10(Mach))*Mach*np.sqrt(T_sg*phy_const.e/m_sg)
                else:
                    psiLarge   = 10**(logPsilarge)
                    sigmaLarge = 10**(f(logPsilarge))
                    functional = np.zeros(np.shape(psiLarge))
                    for k, x in enumerate(Mach*psiLarge):
                        psi_k = psiLarge[k]
                        if x < 1e-2 :
                            series = 0
                            for n in range(6):
                                series += x**(2*n)/((2*n+3)*np.exp(np.sum(np.log(np.linspace(1,2*n+1,2*n+1)))))
                            functional[k] = 2*series*np.exp(-(psi_k**2+Mach**2)/2) #!!! modified !
                        elif x > 1e2:
                            functional[k] = ((1-1/x)*np.exp(-(psi_k-Mach)**2/2)+(1+1/x)*np.exp(-(psi_k+Mach)**2/2))/x**2
                        else:
                            functional[k] = ((1-1/x)*np.exp(psi_k*Mach)+(1+1/x)*np.exp(-psi_k*Mach))*np.exp(-psi_k**2/2-Mach**2/2)/x**2
                        functional[k] *= psi_k**5*sigmaLarge[k]
                    self.K_fixedT[i] = m_sg/self.mass* (phy_const.e*T_sg/(2*np.pi*m_sg))**(1./2.)*integrate.trapz(functional, x= psiLarge)

            print("saving the rates for ", self.name," in ",fileName)
#            with open(fileName,'wb') as f:
#                pkl.dump(self.K_fixedT, f)
            if save:
                print("saving the rates for ", self.name," in ",fileName)
                np.savetxt(fileName, self.K_fixedT)


    def computeK_Momentum(self, u_vec, T_s, T_g): # u is a vector
        m_g       = self.mass + phy_const.m_e
        m_i       = self.mass 
        massRatio = self.mass/m_g
        m_sg      = self.mass/(1 + massRatio)
        T_sg      = (m_g * T_s + m_i * T_g )/(m_i + m_g) 

        psi       = np.sqrt(m_sg*self.Velocity**2/(2*phy_const.e*T_sg))

        Mach_vec      = np.sqrt(m_sg/(2*phy_const.e*T_sg))*u_vec

        logPsilarge    = np.linspace(np.log10(psi[1]), np.log10(psi[-1]), 1000) # We need more resolution of the cross section
        logSigma       = np.log10(self.Sigma)
        logPsi         = np.log10(psi)
        logPsi[0]      = -16
        
        f              = interpolate.interp1d(logPsi, logSigma) #, kind='quadratic')
        result = np.zeros_like(u_vec)
        for i_Mach, Mach in enumerate(Mach_vec):
            u = u_vec[i_Mach]
            if Mach > 30:
                return m_sg/self.mass* 10**f(np.log10(Mach))*u      # Dirac approximation
            else:
                psiLarge   = 10**(logPsilarge)
                sigmaLarge = 10**(f(logPsilarge))
                functional = np.zeros_like(psiLarge)
                for i_y, y in enumerate(2 * Mach * psiLarge):
                    psi_i = psiLarge[i_y]
                    if y < 1e-2 :
                        # Old implementation
                        # series = 0.
                        # for n in range(4):  # Taylor expansion of Eq. (31) in Benilov 97
                        #     series += y**(2*n)/((2*n + 3)*np.exp(np.sum(np.log(np.linspace(1,2*n + 1, 2*n + 1)))))
                        # functional[i_y] = 2*series*np.exp(-(psi_i**2 + Mach**2)) 

                        functional[i_y] = (1./3. + y**2/30 + y**4/840)*np.exp(-(psi_i**2 + Mach**2))*2                  # we multiply to have 16
                        
                    else:          # cosh and sinh give problems
                        functional[i_y] = ((1 + y)*np.exp(-(psi_i + Mach)**2) - (1 - y)*np.exp(-(psi_i - Mach)**2))/(y**3)
                                            
                functional = psiLarge**5*sigmaLarge*functional

            result[i_Mach] = 8*m_sg/self.mass * (phy_const.e*T_sg/(2*np.pi*m_sg))**(1./2.) * integrate.trapz(functional, x= psiLarge)
        
        return result



    def computeK_fct_uT(self, N_T=10, N_u=100, N_integration=1000, N_u_large=10000, N_T_large=1000, sigma_cst=False, save=False):
        m_g = self.mass + phy_const.m_e
        massRatio = self.mass/m_g
        m_sg   = self.mass/(1 + massRatio)
        T_g = 0.025
        if sigma_cst == True:
            fileName = os.path.splitext(self.fileName)[0]+"_K_fct_uT_sigma_cst.txt"
        else:
            fileName = os.path.splitext(self.fileName)[0]+"_K_fct_uT.txt"
        gamma_g = m_g/(phy_const.e*T_g)
        listlogT_s = np.linspace(-2, 1, N_T)
        listT_s   = 10**listlogT_s#+T_g*massRatio)/(1+massRatio)
        #self.Temperature = listT_s
        listU = self.Velocity #np.sqrt(2*sigma[:,0]*phy_const.e/m_sg)
        logU = np.log10(listU)
        logU[0]      = -12
        listlogu_sg = np.linspace(logU[1], logU[-1]-0.5, N_u-1)
        listu_sg = np.array([0]+list(10**listlogu_sg))
        #self.Mach = listu_sg
        list_u_large = np.linspace(0, listu_sg[-1], N_u_large)
        self.Mach = list_u_large
        list_T_large = np.linspace(listT_s[0], listT_s[-1], N_T_large)
        self.Temperature = list_T_large
        if os.path.exists(fileName):
            print("reading the rates for ", self.name," in ",fileName)
    #            with open(fileName,'rb') as file:
    #                self.K_fixedT = pkl.load(file)
            self.K_fct_uT = np.loadtxt(fileName)
        else:
            self.K_fct_uT = np.zeros((N_T, N_u))
            logSigma       = np.log10(self.Sigma)
            f = interpolate.interp1d(logU, logSigma)#, kind='quadratic')
            for i, T_s in enumerate(listT_s):
                gamma_s = self.mass/(phy_const.e*T_s)
                gamma_sg = gamma_s*gamma_g/(gamma_s+gamma_g)
                print('computing K(u,T) for T = %s'%T_s)
                loguLarge = np.linspace(np.log10(listU[1]), np.log10(listU[-1]), N_integration)
                sigmaLarge = 10**(f(loguLarge))
                functional = np.zeros(N_integration)
                psiLarge = np.sqrt(gamma_sg)*10**loguLarge
                for j, u_sg in enumerate(listu_sg):
                    Mach = np.sqrt(gamma_sg)*u_sg
                    if Mach > 30:
                        self.K_fct_uT[i,j] = 0.5*10**f(np.log10(u_sg))*u_sg
                    else:
                        for k, psi_k in enumerate(psiLarge):
                            x = psi_k*Mach
                            if x < 1e-2 :
                                series = 0
                                for n in range(6):
                                    series += x**(2*n)/((2*n+3)*np.exp(np.sum(np.log(np.linspace(1,2*n+1,2*n+1)))))
                                functional[k] = 2*series*np.exp(-(psi_k**2+Mach**2)/2)
                            elif x > 1e2:
                                functional[k] = ((1-1/x)*np.exp(-(psi_k-Mach)**2/2)+(1+1/x)*np.exp(-(psi_k+Mach)**2/2))/x**2
                            else:
                                functional[k] = ((1-1/x)*np.exp(x)+(1+1/x)*np.exp(-x))*np.exp(-psi_k**2/2-Mach**2/2)/x**2
                        self.K_fct_uT[i,j] = 0.5*1/np.sqrt(2*np.pi*gamma_sg)*integrate.trapz(functional*psiLarge**5*sigmaLarge, x=psiLarge)
                #self.K_fct_uT[i,0] = 0.5 * 2/3 * 1/np.sqrt(2*np.pi*gamma_sg) * integrate.trapz(psiLarge**5*np.exp(-psiLarge**2/2)*sigmaLarge, x=psiLarge)
            f_new = interpolate.interp1d(listu_sg, self.K_fct_uT, kind='quadratic')
            self.K_fct_uT = f_new(list_u_large)
            f_new2 = interpolate.interp1d(listT_s, self.K_fct_uT, kind='quadratic', axis=0)
            self.K_fct_uT = f_new2(list_T_large)
            self.K_fct_uT = np.transpose(self.K_fct_uT)
#            with open(fileName,'wb') as f:
#                pkl.dump(self.K_fixedT, f)
            if save:
                print("saving the rates for ", self.name," in ", fileName)
                np.savetxt(fileName, self.K_fct_uT)



    def computeK2_fct_uT(self, N_T=10, N_u=100, N_integration=1000, N_u_large=10000, N_T_large=1000, sigma_cst=False, save=False):
        m_g = self.mass + phy_const.m_e
        massRatio = self.mass/m_g
        m_sg   = self.mass/(1 + massRatio)
        T_g = 0.025
        if sigma_cst == True:
            fileName = os.path.splitext(self.fileName)[0]+"_K2_fct_uT_sigma_cst.txt"
        else:
            fileName = os.path.splitext(self.fileName)[0]+"_K2_fct_uT.txt"
        gamma_g = m_g/(phy_const.e*T_g)
        listlogT_s = np.linspace(np.log10(0.01), 1, N_T)
        listT_s = 10**listlogT_s
        self.Temperature = listT_s
    #    listT_sg   = (10**listlogT_s+T_g*massRatio)/(1+massRatio)
        listU = self.Velocity
        logU = np.log10(listU)
        logU[0] = -12
        listlogu_sg = np.linspace(logU[1], logU[-1]-0.5, N_u-1)
        listu_sg = np.array([0]+list(10**listlogu_sg))
        self.Mach = listu_sg
        list_u_large = np.linspace(0, listu_sg[-1], N_u_large)
        self.Mach = list_u_large
        list_T_large = np.linspace(listT_s[0], listT_s[-1], N_T_large)
        self.Temperature = list_T_large
        if os.path.exists(fileName):
            print("reading the rates for ", self.name," in ",fileName)
    #            with open(fileName,'rb') as file:
    #                self.K_fixedT = pkl.load(file)
            self.K2_fct_uT = np.loadtxt(fileName)
        else:
            self.K2_fct_uT = np.zeros((N_T, N_u))
            logSigma       = np.log10(self.Sigma)
            f = interpolate.interp1d(logU, logSigma)#, kind='quadratic')
            for i, T_s in enumerate(listT_s):
                print('computing K(u,T) for T = %s'%T_s)
                gamma_s = self.mass/(phy_const.e*T_s)
                gamma_sg = gamma_s*gamma_g/(gamma_s+gamma_g)
                loguLarge = np.linspace(np.log10(listU[1]), np.log10(listU[-1]), N_integration)
                sigmaLarge = 10**(f(loguLarge))
                functional = np.zeros(N_integration)
                psiLarge = np.sqrt(gamma_sg)*10**loguLarge
                for j, u_sg in enumerate(listu_sg):
                    Mach = np.sqrt(gamma_sg)*u_sg
                    if Mach > 30:
                        self.K2_fct_uT[i,j] = 0.25 * 10**f(np.log10(u_sg)) * u_sg**3#!!!
                    else:
                        for k, psi_k in enumerate(psiLarge):
                            x = psi_k*Mach
                            if x < 1e-2 :
                                series = 0
                                for n in range(6):
                                    series += x**(2*n)/np.exp(np.sum(np.log(np.linspace(1,2*n+1,2*n+1))))
                                functional[k] = 0.5*series*np.exp(-(psi_k**2+Mach**2)/2)
        #                    elif x > 1e2:
        #                        functional[k] = ((1-1/x)*np.exp(-(psi_k-Mach)**2/2)+(1+1/x)*np.exp(-(psi_k+Mach)**2))/x**2
                            else:
                                functional[k] = 0.25*(np.exp(-(psi_k-Mach)**2/2) - np.exp(-(psi_k+Mach)**2/2))/(psi_k*Mach)
        #                        functional[k] = ((1-1/x)*np.exp(x)+(1+1/x)*np.exp(-x))*np.exp(-psi_k**2/2-Mach**2/2)/x**2
        #                    functional[k] *= psi_k**5*sigmaLarge[k]
                        self.K2_fct_uT[i,j] = 1/np.sqrt(2*np.pi*gamma_sg)*integrate.trapz(functional*psiLarge**5*sigmaLarge, x=psiLarge)*1/ gamma_sg#!!! * 1/u_sg**2
            f_new = interpolate.interp1d(listu_sg, self.K2_fct_uT, kind='quadratic')
            self.K2_fct_uT = f_new(list_u_large)
            f_new2 = interpolate.interp1d(listT_s, self.K2_fct_uT, kind='quadratic', axis=0)
            self.K2_fct_uT = f_new2(list_T_large)
            self.K2_fct_uT = np.transpose(self.K2_fct_uT)
#            with open(fileName,'wb') as f:
#                pkl.dump(self.K2_fct_uT, f)
            if save:
                print("saving the rates for ", self.name," in ", fileName)
                np.savetxt(fileName, self.K2_fct_uT)












































    def computeK_fct_uT_vect(self, N_T=10, N_u=100, N_integration=1000, N_u_large=10000, N_T_large=1000, sigma_cst=False, save=False):
        m_g = self.mass + phy_const.m_e
        massRatio = self.mass/m_g
        m_sg   = self.mass/(1 + massRatio)
        T_g = 0.025
        if sigma_cst == True:
            fileName = os.path.splitext(self.fileName)[0]+"_K_fct_uT_sigma_cst.txt"
        else:
            fileName = os.path.splitext(self.fileName)[0]+"_K_fct_uT.txt"
        gamma_g = m_g/(phy_const.e*T_g)
        listlogT_s = np.linspace(-2, 1, N_T)
        listT_s   = 10**listlogT_s#+T_g*massRatio)/(1+massRatio)
        #self.Temperature = listT_s
        listU = self.Velocity #np.sqrt(2*sigma[:,0]*phy_const.e/m_sg)
        logU = np.log10(listU)
        logU[0]      = -12
        listlogu_sg = np.linspace(logU[1], logU[-1]-0.5, N_u-1)
        listu_sg = np.array([0]+list(10**listlogu_sg))
        #self.Mach = listu_sg
        list_u_large = np.linspace(0, listu_sg[-1], N_u_large)
        self.Mach = list_u_large
        list_T_large = np.linspace(listT_s[0], listT_s[-1], N_T_large)
        self.Temperature = list_T_large
        if os.path.exists(fileName):
            print("reading the rates for ", self.name," in ",fileName)
    #            with open(fileName,'rb') as file:
    #                self.K_fixedT = pkl.load(file)
            self.K_fct_uT = np.loadtxt(fileName)
        else:
            self.K_fct_uT = np.zeros((N_T, N_u))
            logSigma       = np.log10(self.Sigma)
            f = interpolate.interp1d(logU, logSigma)#, kind='quadratic')
            loguLarge = np.linspace(np.log10(listU[1]), np.log10(listU[-1]), N_integration) # velocity list for integration
            sigmaLarge = 10**(f(loguLarge))
            gamma_s = self.mass/(phy_const.e*listT_s)
            gamma_sg = gamma_s*gamma_g/(gamma_s+gamma_g)
            functional = np.zeros(N_T, N_u, N_integration)
            psiLarge = np.sqrt(gamma_sg[:, np.newaxis]).dot(10**loguLarge[np.newaxis, :])
            Mach = np.sqrt(gamma_sg)*listu_sg

            for i, T_s in enumerate(listT_s):
                gamma_s = self.mass/(phy_const.e*T_s)
                gamma_sg = gamma_s*gamma_g/(gamma_s+gamma_g)
                print('computing K(u,T) for T = %s'%T_s)
                loguLarge = np.linspace(np.log10(listU[1]), np.log10(listU[-1]), N_integration)
                sigmaLarge = 10**(f(loguLarge))
                functional = np.zeros(N_integration)
                psiLarge = np.sqrt(gamma_sg)*10**loguLarge
                for j, u_sg in enumerate(listu_sg):
                    Mach = np.sqrt(gamma_sg)*u_sg
                    if Mach > 30:
                        self.K_fct_uT[i,j] = 0.5*10**f(np.log10(u_sg))*u_sg
                    else:
                        for k, psi_k in enumerate(psiLarge):
                            x = psi_k*Mach
                            if x < 1e-2 :
                                series = 0
                                for n in range(6):
                                    series += x**(2*n)/((2*n+3)*np.exp(np.sum(np.log(np.linspace(1,2*n+1,2*n+1)))))
                                functional[k] = 2*series*np.exp(-(psi_k**2+Mach**2)/2)
                            elif x > 1e2:
                                functional[k] = ((1-1/x)*np.exp(-(psi_k-Mach)**2/2)+(1+1/x)*np.exp(-(psi_k+Mach)**2/2))/x**2
                            else:
                                functional[k] = ((1-1/x)*np.exp(x)+(1+1/x)*np.exp(-x))*np.exp(-psi_k**2/2-Mach**2/2)/x**2
                        self.K_fct_uT[i,j] = 0.5*1/np.sqrt(2*np.pi*gamma_sg)*integrate.trapz(functional*psiLarge**5*sigmaLarge, x=psiLarge)
                #self.K_fct_uT[i,0] = 0.5 * 2/3 * 1/np.sqrt(2*np.pi*gamma_sg) * integrate.trapz(psiLarge**5*np.exp(-psiLarge**2/2)*sigmaLarge, x=psiLarge)
            f_new = interpolate.interp1d(listu_sg, self.K_fct_uT, kind='quadratic')
            self.K_fct_uT = f_new(list_u_large)
            f_new2 = interpolate.interp1d(listT_s, self.K_fct_uT, kind='quadratic', axis=0)
            self.K_fct_uT = f_new2(list_T_large)
            self.K_fct_uT = np.transpose(self.K_fct_uT)
#            with open(fileName,'wb') as f:
#                pkl.dump(self.K_fixedT, f)
            if save:
                print("saving the rates for ", self.name," in ", fileName)
                np.savetxt(fileName, self.K_fct_uT)
