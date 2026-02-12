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
#import seaborn as sns
import pandas as pd # type: ignore
from scipy.spatial import Voronoi, voronoi_plot_2d
import json
import os.path
from scipy.optimize import curve_fit
from scipy.spatial import Delaunay
import alphashape # type: ignore
import matplotlib.animation as animation
import celldiv # type: ignore
from scipy.fft import ifft2
from scipy.fft import fft2
from builtins import range, str, int, len, enumerate, min, max, print, open
sys.path.append('../Cell_Shape_Index')
from shape_index import Cell_Shape_Parameter # type: ignore


Files = ['J28_2types_4clusters_set3.xyz']
argv = sys.argv
if argv[1]:
    Files = [argv[1]]

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
file_number=1


def Cell_shape_alignment(f ,th):
    for filename in Files:
        with celldiv.TrajHandle(filename) as th:
            frameCount = 0
            f_count=0
            try:
                for i in range(int(th.maxFrames/nSkip)+1): # i for each frame written to file
                    frameCount += 1

                    if frameCount > args.num_frames:
                        break
                    
                    f=th.ReadFrame(inc=nSkip ) # f is the list of all the cells in the frame
                    f_count+= len(f)
                    m = np.vstack(th.cellInd)
                    f = [ f[i][:180] for i in range(len(m)) ]

                    if len(args.inds) > 0:
                        pass
                    
                    # Nearest neighbor distribution for soft, hard, and all cells
                    m = np.vstack(th.cellInd)
                    all_coords = np.vstack(f)
                    all_coords = np.vstack(all_coords)
                    all_com = []
                    for mi in range(int(len(all_coords)/180)):
                        all_com.append(np.mean(all_coords[mi*180:(mi+1)*180],axis=0))
                    all_com = np.array(all_com)[: , 0:2] #2D coordinates
                    del_tri = Delaunay(all_com)
                    q_array = []
                    area_array = []
                    for coords in all_com[del_tri.simplices]:
                        #shape tensor
                        a = coords[0]
                        b = coords[1]
                        c = coords[2]

                        #check counter clockwise
                        cross_product = (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0])
                        if cross_product < 0:
                            b , c = c , b

                        s1 = [[ b[0]-a[0] , c[0]-a[0]  ]
                            ,[ b[1]-a[1] , c[1]-a[1]  ]]
                        s2 = np.linalg.inv([[1 , 1/2],
                                            [0 , np.sqrt(3)/2]])
                        s = np.dot(s1 , s2)
                        s_tilda = 1/2 * np.array([[ s[0,0] - s[1,1] , s[0,1] + s[1,0] ] , [  s[0,1] + s[1,0] , s[1,1] - s[0,0] ]]) #symmetric traceless matrix
                        theta = np.arctan2(s[1,0] - s[0,1] , s[0,0] + s[1,1])
                        R_neg_theta = [ [ np.cos(theta) , np.sin(theta)] , [ -np.sin(theta) , np.cos(theta)] ]
                        # mag_s_tilda = np.sqrt( np.power(s[0,0],2) + np.power(s[0,1],2) )
                        mag_s_tilda = np.sqrt( s_tilda[0,0]**2 + s_tilda[0,1]**2 )
                        #Division by zero check 
                        if np.linalg.det(s) < 1e-10 or mag_s_tilda < 1e-10:
                            continue
                        q = 1/mag_s_tilda * (np.arcsinh( mag_s_tilda / np.sqrt(np.linalg.det(s)) )) * np.dot(s_tilda , R_neg_theta)
                        q_array.append(q)
                        area = np.linalg.det(s)
                        # print(area)
                        area_array.append(area)

                    # print(q_array)
                    q_area_sum = [ np.array(q_array[i])*area_array[i] for i in range(len(q_array)) ]
                    Q = np.sum(q_area_sum , axis =0) / np.sum(area_array) #check
                    # print(area_array)
                    mag_Q = np.sqrt( np.power( Q[0,0], 2) + np.power( Q[0,1], 2) )
                    Q_s_array = [ np.sqrt( np.power( q_array[i][0,0], 2) + np.power( q_array[i][0,1], 2) ) for i in range(len(q_array)) ]
                    qs_area_sum = [np.array(Q_s_array[i]) * area_array[i] for i in range(len(Q_s_array))] # Is it weighted average?
                    mag_Q_s = np.sum(qs_area_sum , axis =0) / np.sum(area_array) #check
                    Q_alignment = mag_Q/mag_Q_s
                    print( 'Q: '  , mag_Q )
                    print( 'Q_s: ', mag_Q_s)
                    print( 'Q_a: ', mag_Q/mag_Q_s)
                    return Q_alignment


            except celldiv.IncompleteTrajectoryError:
                print ("Stopping...") 


def write_P_Q_toFile(file_name , size , frame_number , f , th):

    Json_file = 'P_Q_{}.json'.format(file_name[:-4])
    Q = Cell_shape_alignment(f,th)
    p_soft, p_hard , p_all = Cell_Shape_Parameter(f , frame_number , file_name , size , plot_system = False , Plot_distribution = False , bin_n=50)
    time = 0 #what was it supposed to be? frame_number*1000?
    if os.path.isfile(Json_file):
        with open(Json_file) as file:
            data = json.load(file)
        data['Q'].append(Q)
        data['p_Hard'].append(p_hard)
        data['p_Soft'].append(p_soft)
        data['p_All'].append(p_all)
        data['Time'].append(time)
        with open(Json_file, 'w') as file:
            json.dump(data, file)
    else:
        data = {
            'Q': [Q],
            'p_Soft': [p_soft],
            'p_Hard': [p_hard],
            'p_All': [p_all],
            'Time': [time]
        }
        with open(Json_file, 'w') as file:
            json.dump(data, file)