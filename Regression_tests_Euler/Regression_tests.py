import sys
from mpi4py import MPI
import muffin                                   ## module in C++ doing the simulation
import glob                                     ## module for unix pathname pattern expansioni
import matplotlib.pyplot as plt                 ## plots
import os                                       ## operating system interface
import numpy as np                              ## Mathematic operations:w
import shutil                                   ## Remove folder
import Euler_tests as Euler

from contextlib import contextmanager
import sys, os

@contextmanager
def suppress_stdout():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:  
            yield
        finally:
            sys.stdout = old_stdout

def bar(t1, t2, tolerance):
    off = abs((np.array(t1)-np.array(t2))).sum()
    return True if off <= tolerance else False

def check_Test(Res, Res_ref, number):

    passed = True
    for i_res, Res in enumerate(Res):
        passed *= bar(Res, Res_ref[i_res], 1e-2)
    if passed ==  True:
        print("TEST "+ str(number)+"\t:\t PASSED")
    else:
        print("TEST "+ str(number)+"\t:\t FAILED")

def plotResult(ax, Directory, color='k', label='test 1'):

    import matplotlib.pyplot as plt

    plt.style.use('classic')
    plt.rcParams["font.family"] = 'Times New Roman'
    plt.rcParams["font.weight"] = 'normal'
    plt.rcParams['figure.facecolor'] = 'white'

    # open all the files in the directory and sort them to do the video in order
    files       = glob.glob(Directory + "/*.txt")
    files.sort(key=os.path.getmtime)
    filename = files[-1]

    results = []
    with open(filename, 'r') as data:
        for line in data:
            p = line.split()
            results.append(np.array(p))

    # Transpose and change data type
    results = np.array(results)
    resultsTP = np.transpose(results)

    resultsTP   = resultsTP.astype(np.float_)        

    # Initializing without the results to choose it directly by hand
    #plotData2 = plotResults.Data()
    def computeEulerVars(resultsArray):
        x       = resultsArray[0]
        rho     = resultsArray[1]
        rhoU    = resultsArray[1]
        rhoE    = resultsArray[2]
        U       = rhoU/rho
        p       = 2./3.*(rhoE - 0.5*rho*U**2)
        T       = p/rho
        return x, rho, rhoU, rhoE, U, p, T

    x, rho, rhoU, rhoE, U, p, T = computeEulerVars(resultsTP)
    ax.plot(x, rho, '--', color=color, linewidth=1.8 , label =label)

# Prepare the plot 
import matplotlib
cmap = matplotlib.cm.get_cmap('plasma')
colors = cmap(np.linspace(1,0,11))

f, ax = plt.subplots(1)
ax.set_xlabel(r'$x$', fontsize=18, weight = 'bold')
ax.set_ylabel(r'$\rho$', fontsize=18)
ax.xaxis.set_tick_params(which='both', size=5, width=1.5, labelsize=13)
ax.yaxis.set_tick_params(which='both', size=5, width=1.5, labelsize=13)
ax.set_xlim([0.3, 0.8])
f.suptitle('Density', fontname = 'Times New Roman',fontsize=16, weight = 'bold')

ax.grid(True)
ax.set_ylim((0, 1.1))


with suppress_stdout():
    ## Test 1
    Res_1_1, Dir = Euler.Euler_first_order("RoeEuler1D")
    Reference_Res_1_1 = [-0.89287,	-0.78000,	-0.51032]
    plotResult(ax, Dir, colors[0],'Test 1.1')
    shutil.rmtree(Dir)

    Res_1_2, Dir = Euler.Euler_first_order("LaxFriedrich")
    Reference_Res_1_2 = [-0.92216,	-0.79202,	-0.51643]
    plotResult(ax, Dir, colors[1],'Test 1.2')
    shutil.rmtree(Dir)

    Res_1_3, Dir = Euler.Euler_first_order("HLLEuler1D")
    Reference_Res_1_3 = [-0.90938,	-0.78780,	-0.51377]
    plotResult(ax, Dir, colors[2],'Test 1.3')
    shutil.rmtree(Dir)

    ## Test 2
    Res_2, Dir = Euler.Euler_TVDRK3()
    Reference_Res_2 = [1.38508, 	1.48389,	1.75175]
    plotResult(ax, Dir, colors[3],'Test 2')
    shutil.rmtree(Dir)


    ## Test 3
    Res_3, Dir = Euler.Euler_SSPRK3()
    Reference_Res_3 = [1.38509, 	1.48389, 1.75175]
    plotResult(ax, Dir, colors[4],'Test 3')
    shutil.rmtree(Dir)

    ## Test 4
    Res_4, Dir = Euler.Euler_SSPRK5()
    Reference_Res_4 = [1.08357,	1.18260,	1.45047]
    plotResult(ax, Dir, colors[5],'Test 4')
    shutil.rmtree(Dir)

    ## Test 5
    Res_5, Dir = Euler.Euler_Roe_second_order('vanAlbada')
    Reference_Res_5 = [ -0.70267,	-0.55467,	-0.28897 ]
    plotResult(ax, Dir, colors[6],'Test 5')
    shutil.rmtree(Dir)

    ## Test 6
    Res_6, Dir = Euler.Euler_Roe_second_order('Venkatakrishnan')
    Reference_Res_6 = [ -0.75499,	-0.59907,	-0.32800 ]
    plotResult(ax, Dir, colors[7],'Test 6')
    shutil.rmtree(Dir)

    ## Test 7
    Res_7, Dir = Euler.Euler_Roe_second_order('ospre')
    Reference_Res_7 = [  -0.68803,	-0.54282,	-0.27834 ]
    plotResult(ax, Dir, colors[8],'Test 7')
    shutil.rmtree(Dir)

    ## Test 8
    Res_8, Dir = Euler.Euler_Roe_second_order('thirdOrder')
    Reference_Res_8 = [ -0.67716,	-0.54125,	-0.27908 ]
    plotResult(ax, Dir, colors[9],'Test 8')
    shutil.rmtree(Dir)

    ## Test 9
    Res_9, Dir = Euler.Euler_first_order_pythonModel()
    Reference_Res_9 = [-0.92216,	-0.79202,	-0.51643]
    plotResult(ax, Dir, colors[10],'Test 9')
    shutil.rmtree(Dir)





check_Test(Res_1_1, Reference_Res_1_1, "1.1")
check_Test(Res_1_2, Reference_Res_1_2, "1.2")
check_Test(Res_1_3, Reference_Res_1_3, "1.3")
check_Test(Res_2, Reference_Res_2, 2)
check_Test(Res_3, Reference_Res_3, 3)
check_Test(Res_4, Reference_Res_4, 4)
check_Test(Res_5, Reference_Res_5, 5)
check_Test(Res_6, Reference_Res_6, 6)
check_Test(Res_7, Reference_Res_7, 7)
check_Test(Res_8, Reference_Res_8, 8)
check_Test(Res_9, Reference_Res_9, 9)

ax.legend()
plt.savefig("./Results.pdf", bbox_inches='tight')

plt.close()


