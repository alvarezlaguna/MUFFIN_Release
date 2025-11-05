#!/usr/bin/env python3
import glob                             ## module for unix pathname pattern expansion
import matplotlib.pyplot as plt         ## plots
import os                               ## operating system interface
import numpy as np                      ## Mathematic operations
import tkinter as tk                     ## To open the file dialog
from tkinter import filedialog

class Data:
    def __init__(self, array=None, filename=None):

        if array is None:
            if filename is None:
                root = tk.Tk()
                root.withdraw()

                currentPath = os.getcwd()
                filename = filedialog.askopenfilename(parent=root, initialdir=currentPath)
                root.destroy()
                print('Reading %s' % filename)

            results = []
            with open(filename, 'r') as data:
                j = 0                   # counter of the lines of the file
                for line in data:
                    p = line.split()
                    results.append(np.array(p))

                # Transpose and change data type
                results = np.array(results)
                resultsTP = np.transpose(results)

                self.resultsArray   = resultsTP.astype(np.float64)        #
                self.nbEqs          = self.resultsArray.shape[0] - 1    # We remove the x axis
                self.nbCells        = self.resultsArray.shape[1]        # We take the number of rows
        else:
            self.resultsArray   = array
            self.nbEqs          = self.resultsArray.shape[0] - 1    # We remove the x axis
            self.nbCells        = self.resultsArray.shape[1]        # We take the number of rows



def getResultsSingleFile(options):
        import tkinter as tk                      ## To open the file dialog
        from tkinter import filedialog

        root = tk.Tk()
        root.withdraw()

        resultsPath = options["resultDir"]
        filename = filedialog.askopenfilename(parent=root,initialdir=resultsPath)
        root.destroy()
        print('Reading %s' % filename)

        with open(filename, 'r') as data:
                NBCELLS = options["nbCells"][0]
                NBEQS   = options["nbEqs"]
                results = np.zeros((NBEQS + 1, NBCELLS))
                j = 0                   # counter of the lines of the file
                for line in data:
                        p = line.split()
                        for i in range(NBEQS + 1):
                                results[i, j] = p[i]
                        j = j + 1       # counter of the lines of the file
                return results

def makePlot(results):
        NBEQS = np.shape(results)[0] - 1        # because we have the
        x = results[0, :]
        f, axarr = plt.subplots(NBEQS, sharex=True, squeeze=False)
        i = 0
        for var_i in range(NBEQS):
                y = results[var_i + 1, :]
                axarr[var_i,0].plot(x, y, linestyle='--', marker='o', color='k', linewidth=1.5, markersize=3)
                axarr[var_i,0].set_xlabel(r'$x$')
                axarr[var_i,0].set_ylabel('var[%i]'%(var_i))
                ymin = min(y)
                ymax = max(y)
                axarr[var_i,0].set_ylim([0.8*ymin,1.2*ymax])
                axarr[var_i,0].grid(True)
                var_i = var_i + 1

        f.suptitle('Simulation')
        plt.show()


def plotSingleResult(options):

        results = getResultsSingleFile(options)
        makePlot(results)

        
def animationResults(varID = 1, fig=None, ax=None, Directory=None):
    from matplotlib import animation, rc
    import sys
    """
        Directory is the name of the directory, by default is null
        fig and ax are the figure and axis that was initialized
    """
    if Directory == None:
        root = tk.Tk()
        root.withdraw()
        currentPath = os.getcwd()
        Directory = filedialog.askdirectory(parent=root, initialdir=currentPath)
   
    # open all the files in the directory and sort them to do the video in order
    files       = glob.glob(Directory + "/*.txt")
    filesSorted = sorted(files, key = lambda x: os.path.getmtime(x), reverse=True)
    files.sort(key=os.path.getmtime)
    
    if fig == None:
        fig, ax = plt.subplots()

        ax.set_xlim(( 0, 1))
        ax.set_ylim((0, 1.2))
        ax.grid(True)
    line, = ax.plot([], [], lw=2, color='k')
        
    def animate_func(i):
        file_name = files[i]
        # Take the time
        base = base=os.path.basename(file_name)
        filenameBase = os.path.splitext(base)[0]
        time = filenameBase.rsplit('_', 1)[1]
        # Initializing data with the file name
        resultsData = Data(array=None,filename = file_name)
        # Check prints
        if varID > np.shape(resultsData.resultsArray)[0] - 1:
            print("varID not available")
            sys.exit()
        if varID == 0:
            print("varID == 0 is the x variable")
            sys.exit()
        
        x       = resultsData.resultsArray[0, :]
        y       = resultsData.resultsArray[varID, :]
        line.set_data(x, y)
        return line,

    anim = animation.FuncAnimation(
                                    fig,
                                    animate_func,
                                    frames = int(np.shape(files)[0])
                                    )
    plt.close()
    
    return anim
    



