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

            # detect file type by extension and read accordingly
            _, ext = os.path.splitext(filename)
            ext = ext.lower()
            if ext in ('.h5', '.hdf5'):
                try:
                    import h5py
                except Exception as e:
                    raise RuntimeError('h5py required to read HDF5 files: %s' % e)
                with h5py.File(filename, 'r') as f:
                    # Expect datasets 'x' and 'u' (u: NBCELLS x NBEQS or NBEQS x NBCELLS)
                    if 'x' in f:
                        x = np.array(f['x'])
                    else:
                        raise RuntimeError('HDF5 file does not contain "x" dataset')
                    if 'u' in f:
                        u = np.array(f['u'])
                    else:
                        # try alternatives (Phi, dataset per variable) - not supported
                        raise RuntimeError('HDF5 file does not contain "u" dataset')

                    # normalize shapes: we want resultsArray with shape (NBEQS+1, NBCELLS)
                    # if u is shape (NBCELLS, NBEQS) -> transpose to (NBEQS, NBCELLS)
                    if u.ndim != 2:
                        raise RuntimeError('u dataset must be 2D')
                    if u.shape[0] == x.shape[0]:
                        u2 = u.T
                    elif u.shape[1] == x.shape[0]:
                        u2 = u
                    else:
                        # ambiguous orientation: try to flatten/reshape
                        raise RuntimeError('u dataset shape %s is incompatible with x length %d' % (str(u.shape), x.shape[0]))

                    # set time if present
                    if 'time' in f:
                        try:
                            self.time = float(np.array(f['time']))
                        except Exception:
                            # time could be stored as numpy scalar or dataset
                            self.time = float(np.array(f['time']).tolist())
                    else:
                        # try to infer time from filename (last underscore token)
                        try:
                            base = os.path.basename(filename)
                            filenameBase = os.path.splitext(base)[0]
                            self.time = float(filenameBase.rsplit('_', 1)[-1])
                        except Exception:
                            self.time = None

                    # If Phi dataset exists, append it as an extra variable row
                    if 'Phi' in f:
                        phi = np.array(f['Phi'])
                        if phi.size != x.size:
                            raise RuntimeError('Phi dataset length %d incompatible with x length %d' % (phi.size, x.size))
                        self.resultsArray = np.vstack((x, u2, phi))
                    else:
                        self.resultsArray = np.vstack((x, u2))
                    self.nbEqs = self.resultsArray.shape[0] - 1
                    self.nbCells = self.resultsArray.shape[1]
            else:
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
                    # Try to infer time from filename for text files
                    try:
                        base = os.path.basename(filename)
                        filenameBase = os.path.splitext(base)[0]
                        self.time = float(filenameBase.rsplit('_', 1)[-1])
                    except Exception:
                        self.time = None
        else:
            self.resultsArray   = array
            self.nbEqs          = self.resultsArray.shape[0] - 1    # We remove the x axis
            self.nbCells        = self.resultsArray.shape[1]        # We take the number of rows
            self.time = None



def getResultsSingleFile(options):
        import tkinter as tk                      ## To open the file dialog
        from tkinter import filedialog

        root = tk.Tk()
        root.withdraw()

        resultsPath = options["resultDir"]
        filename = filedialog.askopenfilename(parent=root,initialdir=resultsPath)
        root.destroy()
        print('Reading %s' % filename)

        # Support both text and HDF5 files
        _, ext = os.path.splitext(filename)
        ext = ext.lower()
        if ext in ('.h5', '.hdf5'):
            try:
                import h5py
            except Exception as e:
                raise RuntimeError('h5py required to read HDF5 files: %s' % e)
            with h5py.File(filename, 'r') as f:
                if 'x' not in f or 'u' not in f:
                    raise RuntimeError('HDF5 file missing required datasets (x, u)')
                x = np.array(f['x'])
                u = np.array(f['u'])
                # orient u to shape (NBCELLS, NBEQS)
                if u.ndim != 2:
                    raise RuntimeError('u dataset must be 2D')
                if u.shape[0] == x.size:
                    u2 = u
                elif u.shape[1] == x.size:
                    u2 = u.T
                else:
                    raise RuntimeError('u dataset shape %s incompatible with x length %d' % (str(u.shape), x.size))
                NBEQS = u2.shape[1]
                NBCELLS = x.size
                # Build results with same shape as text reader: (NBEQS+1, NBCELLS)
                base = np.vstack((x, u2.T))
                # Append Phi if present
                if 'Phi' in f:
                    phi = np.array(f['Phi'])
                    if phi.size != x.size:
                        raise RuntimeError('Phi dataset length %d incompatible with x length %d' % (phi.size, x.size))
                    results = np.vstack((base, phi)).astype(np.float64)
                else:
                    results = base.astype(np.float64)
                return results
        else:
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
    



