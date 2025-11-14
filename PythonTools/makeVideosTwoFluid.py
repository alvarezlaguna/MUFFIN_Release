import glob                             ## module for unix pathname pattern expansion
import matplotlib
import matplotlib.pyplot as plt         ## plots
import sys
import os                               ## operating system interface
import numpy as np                      ## Mathematic operations

# Prefer importing from installed package (when muffin is installed).
try:
    from PythonTools import plotResults
except Exception:
    # Fallback for local development
    # keep previous relative path behavior only when package not installed
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
    import plotResults


def make_videos(ResultDir: str, MASSRATIO: float = 100000.0, saverate: int = 2, frame_rate: int = 20):
    """Generate PNG frames and MP4 videos from text result files.

    Parameters
    ----------
    ResultDir : str
        Path to directory containing `*.txt` result files.
    MASSRATIO : float
        Mass ratio used to compute reference potentials for plotting.
    saverate : int
        Only process every `saverate`-th file (to reduce frames).
    frame_rate : int
        Frame rate used when invoking ffmpeg to create mp4 files.

    Returns
    -------
    dict
        Paths to the generated mp4 files.
    """
    FigureDir = os.path.join(ResultDir, "Figures")
    os.makedirs(FigureDir, exist_ok=True)

    # collect both text and HDF5 files
    files = glob.glob(os.path.join(ResultDir, "*.txt"))
    files += glob.glob(os.path.join(ResultDir, "*.h5"))
    files += glob.glob(os.path.join(ResultDir, "*.hdf5"))
    files = sorted(files, key=os.path.getmtime)

    frame_idx = 0
    file_idx = 0
    for filename in files:
        base = os.path.basename(filename)
        filenameBase = os.path.splitext(base)[0]
        # assume time is the last underscore-separated token
        time = filenameBase.rsplit('_', 1)[-1]
        file_idx += 1
        if file_idx % saverate == 0:
            frame_idx += 1
            # Initialize plot data from file (supports .txt and .h5/.hdf5)
            plotData = plotResults.Data(array=None, filename=filename)

            plt.rcParams["font.family"] = 'Times New Roman'

            f, ax = plt.subplots(1)

            # Density
            x = plotData.resultsArray[0]
            n_e = plotData.resultsArray[1]
            n_i = plotData.resultsArray[3]

            ax.plot(x, n_i, color=(255/255, 97/255, 3/255), linewidth=1.8, markersize=3, label='Ions')
            ax.plot(x, n_e, linestyle='--', color='k', linewidth=1.8, markersize=3, label='Electrons')
            ax.set_xlabel(r'$x/\lambda$', fontsize=18, weight='bold')
            ax.set_ylabel(r'$n$', fontsize=18)
            ax.xaxis.set_tick_params(which='both', size=5, width=1.5, labelsize=13)
            ax.yaxis.set_tick_params(which='both', size=5, width=1.5, labelsize=13)
            plt.legend(loc='upper left', fontsize=12)

            f.suptitle('Density', fontname='Times New Roman', fontsize=16, weight='bold')
            # prefer time from HDF5 dataset if available
            use_time = getattr(plotData, 'time', None)
            if use_time is None:
                # fallback to filename-based time
                try:
                    base = os.path.basename(filename)
                    filenameBase = os.path.splitext(base)[0]
                    use_time = float(filenameBase.rsplit('_', 1)[-1])
                except Exception:
                    use_time = None
            ax.set_title(r'$t = $' + (str(round(float(use_time), 2)) if use_time is not None else 'N/A'), loc='right')

            # ax.set_ylim([0,2.2])
            ax.set_ylim([0, 2.2])
            ax.grid(True)
            plt.savefig(os.path.join(FigureDir, "density_%02d.png" % frame_idx), dpi=150)
            plt.close(f)

            f, ax = plt.subplots(1)
            nU_e = plotData.resultsArray[2]
            nU_i = plotData.resultsArray[4]
            ax.plot(x, nU_i, color=(255/255, 97/255, 3/255), linewidth=1.8, markersize=3, label='Ions')
            ax.plot(x, nU_e, linestyle='--', color='k', linewidth=1.8, markersize=3, label='Electrons')
            ax.set_xlabel(r'$x/\lambda$', fontsize=18, weight='bold')
            ax.set_ylabel(r'$nu$', fontsize=18)
            ax.xaxis.set_tick_params(which='both', size=5, width=1.5, labelsize=13)
            ax.yaxis.set_tick_params(which='both', size=5, width=1.5, labelsize=13)
            plt.legend(loc='upper left', fontsize=12)

            f.suptitle('Flux', fontname='Times New Roman', fontsize=16, weight='bold')
            ax.set_title(r'$t = $' + (str(round(float(use_time), 2)) if use_time is not None else 'N/A'), loc='right')

            ymin = min(nU_e)
            ymax = max(nU_e)
            ax.set_ylim([1.4 * ymin, 1.4 * ymax])
            ax.grid(True)
            plt.savefig(os.path.join(FigureDir, "flux_%02d.png" % frame_idx), dpi=150)
            plt.close(f)

            f, ax = plt.subplots(1)
            U_e = nU_e[:] / n_e[:]
            U_i = nU_i[:] / n_i[:]
            ax.plot(x, U_i, color=(255/255, 97/255, 3/255), linewidth=1.8, markersize=3, label='Ions')
            ax.plot(x, U_e, linestyle='--', color='k', linewidth=1.8, markersize=3, label='Electrons')
            ax.set_xlabel(r'$x/\lambda$', fontsize=18, weight='bold')
            ax.set_ylabel(r'$u$', fontsize=18)
            ax.xaxis.set_tick_params(which='both', size=5, width=1.5, labelsize=13)
            ax.yaxis.set_tick_params(which='both', size=5, width=1.5, labelsize=13)
            plt.legend(loc='upper left', fontsize=12)

            f.suptitle('Velocity', fontname='Times New Roman', fontsize=16, weight='bold')
            ax.set_title(r'$t = $' + (str(round(float(use_time), 2)) if use_time is not None else 'N/A'), loc='right')

            ax.set_ylim([-4, 4])
            ax.grid(True)
            plt.savefig(os.path.join(FigureDir, "velocity_%02d.png" % frame_idx), dpi=150)
            plt.close(f)

            f, ax = plt.subplots(1)
            phi = plotData.resultsArray[5]

            phi_W = np.log((MASSRATIO / (2 * np.pi)) ** (1 / 2)) * np.ones(x.size)
            phi_P = (phi_W + 0.5) * np.ones(x.size)

            ax.plot(x, phi_W, color=(0, 0.4, 0), linewidth=1.8, markersize=3)
            ax.plot(x, phi_P, color=(0, 0.6, 0), linewidth=1.8, markersize=3)
            ax.plot(x, phi, color=(61/255, 89/255, 171/255), linewidth=1.8, markersize=3, label='Potential')

            ax.set_xlabel(r'$x/\lambda$', fontsize=18, weight='bold')
            ax.set_ylabel(r'$\phi$', fontsize=18)
            ax.xaxis.set_tick_params(which='both', size=5, width=1.5, labelsize=13)
            ax.yaxis.set_tick_params(which='both', size=5, width=1.5, labelsize=13)

            f.suptitle('Potential', fontname='Times New Roman', fontsize=16, weight='bold')
            ax.set_title(r'$t = $' + (str(round(float(use_time), 2)) if use_time is not None else 'N/A'), loc='right')

            plt.text(0.9 * max(x), 1.04 * max(phi_W), r'$\phi_W$', fontsize=13, color=(0, 0.4, 0), weight='bold')
            plt.text(0.9 * max(x), 1.04 * max(phi_P), r'$\phi_p$', fontsize=13, color=(0, 0.6, 0), weight='bold')

            ax.set_ylim([0, 7])
            ax.grid(True)
            plt.savefig(os.path.join(FigureDir, "potential_%02d.png" % frame_idx), dpi=150)
            plt.close(f)

    # Build videos with ffmpeg
    density_mp4 = os.path.join(ResultDir, "density.mp4")
    flux_mp4 = os.path.join(ResultDir, "flux.mp4")
    velocity_mp4 = os.path.join(ResultDir, "velocity.mp4")
    potential_mp4 = os.path.join(ResultDir, "potential.mp4")

    os.system(f"ffmpeg -r {frame_rate} -i {os.path.join(FigureDir, 'density_%02d.png')} -vcodec mpeg4 -y {density_mp4}")
    os.system(f"ffmpeg -r {frame_rate} -i {os.path.join(FigureDir, 'flux_%02d.png')} -vcodec mpeg4 -y {flux_mp4}")
    os.system(f"ffmpeg -r {frame_rate} -i {os.path.join(FigureDir, 'velocity_%02d.png')} -vcodec mpeg4 -y {velocity_mp4}")
    os.system(f"ffmpeg -r {frame_rate} -i {os.path.join(FigureDir, 'potential_%02d.png')} -vcodec mpeg4 -y {potential_mp4}")

    return {
        'density': density_mp4,
        'flux': flux_mp4,
        'velocity': velocity_mp4,
        'potential': potential_mp4,
    }


if __name__ == '__main__':
    # Allow calling from command-line: provide folder and MASSRATIO
    import argparse

    parser = argparse.ArgumentParser(description='Create videos from MUFFIN two-fluid results')
    parser.add_argument('folder', nargs='?', default='../testcases/test_Sheath', help='Result folder containing .txt files')
    parser.add_argument('--massratio', type=float, default=100000.0, help='Mass ratio used for potential plotting')
    parser.add_argument('--saverate', type=int, default=2, help='Process every saverate-th file')
    parser.add_argument('--framerate', type=int, default=20, help='Frame rate for output videos')

    args = parser.parse_args()
    make_videos(args.folder, MASSRATIO=args.massratio, saverate=args.saverate, frame_rate=args.framerate)
