import csv
import sys
import os
import argparse
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import multivariate_normal, gaussian_kde
from scipy.spatial import Voronoi, voronoi_plot_2d
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib.colors as mcol
from scipy.spatial import ConvexHull, convex_hull_plot_2d
import pandas as pd
from scipy.spatial import Voronoi, voronoi_plot_2d
import json
import os.path
from scipy.optimize import curve_fit
from scipy.spatial import Delaunay
import alphashape
import matplotlib.animation as animation
import celldiv # type: ignore
from scipy.fft import ifft2
from scipy.fft import fft2

Files = ['J28_2types_4clusters_set3.xyz'] # example file
argv = sys.argv
if argv[1]:
    Files = [argv[1]]


# System units:
# unit length found by using the calculate z function
Unit_length = 1.53/2  # Radius_of_a_cell    


parser = argparse.ArgumentParser()

parser.add_argument("trajPath", type=str,
                    help="Trajectory path. Absolute or relative.")

parser.add_argument("-k", "--skip", type=int, required=False,
                    help="Trajectory frame skip rate. E.g. SKIP=10 will only \
                    render every 10th frame.",
                    default=1)

parser.add_argument("--min-cells", type=int, required=False,
                    help='Start rendering when system has at least this many cells',
                    default=1)

parser.add_argument("--inds", type=int, required=False, nargs='+',
                    help="Only render cells with these indices",
                    default=[])

parser.add_argument("-nf", "--num-frames", type=int, required=False,
                    help="Only render this many frames.",
                    default=sys.maxsize)

args = parser.parse_args(argv[1:])

nSkip = args.skip
# fig = plt.figure()
# fig.set_size_inches(24,20)
# ax = fig.add_subplot()
# ax.set_xlabel('X')
# ax.set_ylabel('Y')

file_number=1
cell_radius = 0.57


from matplotlib import rcParams, rc, ticker, colors, cm
rc('font', **{'family': 'serif', 'serif': ['Times New Roman']})
rc('text', usetex=True)


firstfaces = []
f_count = 0
DP = []
DC = []
time = []


z_contact = 21 #4.5

z_precursor = 19#1.5

# Animate_cell_elongation(view= "XZ" , size_z = 50, frame_rate=5 , z_contact = 0 , z_precursor = 0)