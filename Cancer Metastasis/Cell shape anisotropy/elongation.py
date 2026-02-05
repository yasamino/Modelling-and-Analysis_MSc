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
from builtins import range, str, int, len, enumerate, min, max, print


Files = ['J28_2types_4clusters_set3.xyz']
argv = sys.argv
if argv[1]:
    Files = [argv[1]]

Unit_length = 1.53/2  # Radius_of_a_cell

# mpl.use('pdf')
# plt.rc('font', family='serif', serif='Times')
# plt.rc('text', usetex=True)
# plt.rc('xtick', labelsize=8)
# plt.rc('ytick', labelsize=8)
# plt.rc('axes', labelsize=8)


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

def plot_distribution_of_elongation(array , bins, j , Filename):
    x_max = 1
    x=np.linspace(min(array), max(array), bins+1)
    y=[0]*(bins+1)
    for i in array:
        y[int((i-min(array))/(max(array)-min(array))*bins)]+=1
    plt.scatter(x,np.array(y)/len(array) *100 , label='Data' )
    # gauss = lambda x, sigma2_times_2_rev , mu: np.sqrt(sigma2_times_2_rev/np.pi) * np.exp(-1*sigma2_times_2_rev * (x-mu)**2)
    # exponential = lambda x, a, b: a*np.exp(-b*x)

    # percent = 0.30
    # gauss_x = x[:int(percent*bins)]
    # gauss_y = np.array(y)[:int(percent*bins)]/len(array) *100

    # exp_x = x[int(percent*bins):]
    # exp_y = np.array(y)[int(percent*bins):]/len(array) *100

    # parameter_gauss, covarience_gauss = curve_fit(gauss, gauss_x, gauss_y , maxfev=5000)
    # parameter_exp , covarience_exp = curve_fit(exponential, exp_x, exp_y , maxfev=5000)
    #parameter_gauss_all , covarience_gauss_all = curve_fit(gauss, x, np.array(y)/len(array) *100 , maxfev=5000)

    # y_data_fit_gauss = gauss(gauss_x, parameter_gauss[0], parameter_gauss[1])
    # y_data_fit_exp = exponential(exp_x, parameter_exp[0], parameter_exp[1])
    #y_data_fit_gauss_all = gauss(x, parameter_gauss_all[0], parameter_gauss_all[1])

    #plt.plot(gauss_x, y_data_fit_gauss, label='$\\frac{1}{\sqrt{2 \pi  %5.3f}} e^{- \\frac{(\\Upsilon - %5.3f )^2}{%5.3f}} $' % (1/(parameter_gauss[0]*2) ,parameter_gauss[1] ,1/(parameter_gauss[0]*2)) , color = 'g')
    # plt.plot(gauss_x, y_data_fit_gauss, label='$ \\sigma^2 =  %5.3f , \\mu = %5.3f $'%( 1/(parameter_gauss[0]*2) , parameter_gauss[1]) , color = 'g')
    # plt.plot(exp_x, y_data_fit_exp, label='$%5.3f e^{- %5.3f \\Upsilon} $' % (parameter_exp[0], parameter_exp[1]) , color = 'orange')
    #plt.plot(x, y_data_fit_gauss_all, label='$ fit: \\frac{1}{\sqrt{2 \pi  %5.3f}} e^{- \\frac{(x - %5.3f )^2}{%5.3f}} $' % (1/(parameter_gauss_all[0]*2) ,parameter_gauss_all[1] ,1/(parameter_gauss_all[0]*2)) , color = 'k')

    plt.xlabel('$\\Upsilon$')
    plt.ylabel('N (%)')
    plt.yscale('log')
    plt.legend()
    plt.ylim(0.01,15)
    plt.xlim(0,x_max)
    plt.title('Distribution of $\\Upsilon$ , time = {}'.format(j*10000))
    plt.savefig('Dist_elong_{}_{}.png'.format(Filename, j))
    plt.clf()


def Elongation( Plot_distribution = True , bin_n=50 ):
    for filename in Files:
        f_count=0
        with celldiv.TrajHandle(filename) as th:
            try:
                for frame in range(int(th.maxFrames/nSkip)+1): # i for each frame written to file
                    Elongation_of_soft_cells = []
                    Elongation_of_hard_cells = []
                    if frame+1 > args.num_frames:
                        break
                    
                    f=th.ReadFrame(inc=nSkip ) # f is the list of all the cells in the frame
                    f_count+= len(f)
                    m = np.vstack(th.cellInd)
                    s1 = np.where(m < 0)[0]
                    s2 = np.where(m >= 0 )[0]
                    soft = [ f[i][:180] for i in s1 ]
                    hard = [ f[i][:180] for i in s2 ]
                    
                    for cell in soft:
                        hull = ConvexHull(cell[:,0:2]) # Finding the boundary points
                        Area = hull.volume
                        p = cell[hull.vertices,0:2]
                        COM = np.mean(p, axis=0)
                        centered_p = p - COM

                        # calculations based on paper by Comelles et.al: https://doi.org/10.7554/eLife.57730
                        ''' Exx = D_alpha,  Eyy = -D_alpha,  Exy = Eyx = B_alpha '''
                        '''D_alpha = 1/A * sum( Cos(2theta) dA )'''
                        '''B_alpha = 1/A * sum( Sin(2theta) dA )'''
                        D_alpha = np.float64(0)
                        B_alpha = np.float64(0)
                        Theta = np.arctan2(centered_p[:,1], centered_p[:,0]) # angle of each point with respect to COM
                        for i in range(len(centered_p[:,0])):
                            d_A = 0.5 * np.abs(centered_p[i-1][0] * centered_p[i][1] - centered_p[i][0] * centered_p[i-1][1])
                            D_alpha += np.cos(2*Theta[i]) * d_A
                            B_alpha += np.sin(2*Theta[i]) * d_A
                        D_alpha = D_alpha/Area
                        B_alpha = B_alpha/Area
                        Mag_of_elongation = np.sqrt(D_alpha**2 + B_alpha**2)
                        Elongation_of_soft_cells.append(Mag_of_elongation)
                    
                    for cell in hard:
                        hull = ConvexHull(cell[:,0:2]) # Finding the boundary points
                        Area = hull.volume
                        p = cell[hull.vertices,0:2]
                        COM = np.mean(p, axis=0)
                        centered_p = p - COM

                        # calculations based on paper by Comelles et.al: https://doi.org/10.7554/eLife.57730
                        ''' Exx = D_alpha,  Eyy = -D_alpha,  Exy = Eyx = B_alpha '''
                        '''D_alpha = 1/A * sum( Cos(2theta) dA )'''
                        '''B_alpha = 1/A * sum( Sin(2theta) dA )'''
                        D_alpha = np.float64(0)
                        B_alpha = np.float64(0)
                        Theta = np.arctan2(centered_p[:,1], centered_p[:,0]) # angle of each point with respect to COM
                        for i in range(len(centered_p[:,0])):
                            d_A = 0.5 * np.abs(centered_p[i-1][0] * centered_p[i][1] - centered_p[i][0] * centered_p[i-1][1])
                            D_alpha += np.cos(2*Theta[i]) * d_A
                            B_alpha += np.sin(2*Theta[i]) * d_A
                        D_alpha = D_alpha/Area
                        B_alpha = B_alpha/Area
                        Mag_of_elongation = np.sqrt(D_alpha**2 + B_alpha**2)
                        Elongation_of_hard_cells.append(Mag_of_elongation)
                    
                    
                    if Plot_distribution:
                        plt.clf()
                        plt.hist(Elongation_of_hard_cells, bins=bin_n, density=True, alpha=0.5, label='Hard cells', color='r' )
                        plt.hist(Elongation_of_soft_cells, bins=bin_n, density=True, alpha=0.5, label='Soft cells', color='b' )
                        # x_max = 0.6
                        # print(max(Elongation_of_hard_cells))
                        # x=np.linspace(min(Elongation_of_hard_cells), max(Elongation_of_hard_cells), bin_n+1)
                        # y=[0]*(bin_n+1)
                        # for i in Elongation_of_hard_cells:
                        #     y[int((i-min(Elongation_of_hard_cells))/(max(Elongation_of_hard_cells)-min(Elongation_of_hard_cells))*bin_n)]+=1
                        # plt.scatter(x,np.array(y)/len(Elongation_of_hard_cells) *100 , label='Hard cells', color='r' )
                        # print(max(Elongation_of_soft_cells))
                        # x= np.linspace(min(Elongation_of_soft_cells), max(Elongation_of_soft_cells), bin_n+1)
                        # y=[0]*(bin_n+1)
                        # for i in Elongation_of_soft_cells:
                        #     y[int((i-min(Elongation_of_soft_cells))/(max(Elongation_of_soft_cells)-min(Elongation_of_soft_cells))*bin_n)]+=1
                        # plt.scatter(x,np.array(y)/len(Elongation_of_soft_cells) *100 , label='Soft cells', color='b' )

                        plt.xlabel('$\\Upsilon$')
                        plt.ylabel('N (%)')
                        plt.xlim(0,0.6)
                        plt.savefig('Dist_elong_{}_{}.png'.format(filename[:-4], frame))
                    
            except celldiv.IncompleteTrajectoryError:
                print ("Stopping...") 


def Animate_cell_elongation(view = 'XY', size_z = 28, frame_rate = 5 , z_contact = 3.5 , z_precursor = 1.5):
    start_from =0
    end = 0
    size=100
    array_of_coms_x = []
    array_of_coms_y = []
    #Temporary 
    draw_lines = 0
    # z_contact = 3.5
    # z_precursor = 1.5

    def soft_hard_animation(f,th,frame_number , Filename , size , ax , view = 'XY' , z = 1.9):
        m = np.vstack(th.cellInd)
        c = [ f[i][:180] for i in range(len(m)) ]
        Elongation_mag = []
        index_of_elongated_cells = []
        
        #Elong
        
        for index , cell in enumerate(c):
            # print("cell: ", cell)
            hull = ConvexHull(cell[:,0:2])
            Area = hull.volume
            p = cell[hull.vertices,0:2]
            COM = np.mean(p, axis=0)
            array_of_coms_x.append(COM[0])
            array_of_coms_y.append(COM[1])
            centered_p = p - COM
            '''D_alpha = 1/A * sum( Cos(2theta) dA )'''
            '''B_alpha = 1/A * sum( Sin(2theta) dA )'''
            D_alpha = np.float64(0)
            B_alpha = np.float64(0)
            Theta = np.arctan2(centered_p[:,1], centered_p[:,0])
            for i in range(len(Theta)):
               d_A = 0.5 * np.abs(centered_p[i-1][0] * centered_p[i][1] - centered_p[i][0] * centered_p[i-1][1])
               D_alpha += np.cos(2*Theta[i]) * d_A
               B_alpha += np.sin(2*Theta[i]) * d_A
            D_alpha = D_alpha/Area
            B_alpha = B_alpha/Area
            Magnitude_of_elongation = np.linalg.norm([[D_alpha, B_alpha], [B_alpha, -D_alpha]]) # Norm of the matrix is the magnitude of elongation
            Elongation_mag.append(Magnitude_of_elongation)
            Mag_limit = 0.5
            
            if view == 'XY':
                pass
            elif view == 'YZ':
                hull = ConvexHull(cell[:,1:3])
                p = cell[hull.vertices,1:3]
            elif view == 'XZ':
                hull = ConvexHull(cell[:,::2])
                p = cell[hull.vertices,::2]
            

            patches = []
            patches.append(Polygon(p, closed=True))
            p = PatchCollection(patches,edgecolor="r",alpha=0.2)
            # set all to a nice purple color
            # p = PatchCollection(patches, color=(0.5, 0.2, 0.7, 1.0), edgecolor=(0.5, 0.2, 0.7, 1.0), alpha=0.3)

            p.set_color(cpick.to_rgba(Magnitude_of_elongation))
            ax.add_collection(p)
        # plt.plot(array_of_coms_x , array_of_coms_y)

        if view == 'XY':
            ax.set_xlim([0,size])
            ax.set_ylim([0,size])
        elif view == 'YZ':
            ax.set_xlim([0,size])
            ax.set_ylim([0,z])
        elif view == 'XZ':
            ax.set_xlim([0,size])
            ax.set_ylim([0,z])

        # ax.set_xlim([10,size-10])
        # ax.set_ylim([10,size-10])
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_linewidth(2)

        if draw_lines:
            #draw a horizontal line at z_contact
            ax.axhline(y= z_contact, color='k', linestyle='--', linewidth=1 , label='Contact layer')
            #draw a horizontal line at z_precursor
            ax.axhline(y= z_precursor, color='k', linestyle=':', linewidth=1 , label='Precursor')
            ax.legend(loc='upper left', fontsize=10, frameon=False)
            
        # ax.set_title("Cell Elongation, Time = {}".format(frame_number*100))
        print(frame_number)
        #save fig
        # if frame_number >=0:
        #     plt.savefig('elongation_{}_{}.png'.format(Filename, frame_number))
        return ax.plot()
    
    def update(frame_number):
        ax.clear()
        frame = th.ReadFrame(inc=nSkip)
        return soft_hard_animation(frame, th, frame_number , filename, size, ax , view , size_z) #### set the size

    for filename in Files:
        with celldiv.TrajHandle(filename) as th:
            if view == 'XY':
                fig, ax = plt.subplots(figsize = (size/5, size/5))
            else: # view == 'YZ':
                fig, ax = plt.subplots(figsize=(size/5,size_z/5)) ####### set the figsize for ephitelial tissue

            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_aspect('equal', adjustable='box')
            # range of x and y axis
            ax.set_xlim([25,size-25])
            ax.set_ylim([25,size-25])
            # remove the numbers from the axis

            cm1 = mcol.LinearSegmentedColormap.from_list("MyCmapName",['#2A9B2A','r'])
            cnorm = mcol.Normalize(vmin=0,vmax=1)
            cpick = cm.ScalarMappable(norm=cnorm,cmap=cm1)
            cpick.set_array([])
            #create an animation
            for _ in range(start_from):
                th.ReadFrame(inc=nSkip)
            ani = animation.FuncAnimation(fig=fig , func=update, frames=range(start_from,th.maxFrames-end,nSkip), interval=50, repeat=False)
            cbar = plt.colorbar(cpick, ax=ax , pad =0.1)
            cbar.set_label("$\\Upsilon$", labelpad=10)  # Main label (on the side)
            cbar.ax.set_title('elongated', pad=10, fontsize=12)  # Label above the color bar
            cbar.ax.set_xlabel('Spherical', labelpad=10 , fontsize=12)  # Label below the color bar
            # plt.show()
            #save animation as gif
            if view == 'XY':
                ani.save('elong_colorma_XY_{}.gif'.format(filename[:-4]), writer='Pillow', fps=frame_rate)
            elif view == 'YZ':
                ani.save('elong_colormap_YZ_{}.gif'.format(filename[:-4]), writer='Pillow', fps=frame_rate)
            elif view == 'XZ':
                ani.save('elong_colormap_XZ_{}.gif'.format(filename[:-4]), writer='Pillow', fps=frame_rate)
            print("The view is: {}".format(view))


def Animate_layers_elongation(): #interactive - change the parameters
    start_from = 0
    end = 0
    
    Z = float(input("Size of the box - Z axis : \n"))
    layers = int(input("Number of layers : \n"))
    size = int(input("Size of the box - X and Y : \n"))
    view = input("view : \n")

    def soft_hard_animation(f,th,frame_number , Filename , size , ax , view = 'XY' , z = 1.9):
        m = np.vstack(th.cellInd)
        c = [ f[i][:180] for i in range(len(m)) ]
        Elongation_mag = []
        index_of_elongated_cells = []
        #Elong
        
        for index , cell in enumerate(c):

            COM_Z = np.mean(cell[:,2], axis=0) 
            
            hull = ConvexHull(cell[:,0:2])
            Area = hull.volume
            p = cell[hull.vertices,0:2]
            COM = np.mean(p, axis=0)# # The center of mass is calculated using ONLY the outer point of the cell
            centered_p = p - COM
            '''D_alpha = 1/A * sum( Cos(2theta) dA )'''
            '''B_alpha = 1/A * sum( Sin(2theta) dA )'''
            D_alpha = np.float64(0)
            B_alpha = np.float64(0)
            Theta = np.arctan2(centered_p[:,1], centered_p[:,0])
            for i in range(len(Theta)):
                d_A = 0.5 * np.abs(centered_p[i-1][0] * centered_p[i][1] - centered_p[i][0] * centered_p[i-1][1])
                D_alpha += np.cos(2*Theta[i]) * d_A
                B_alpha += np.sin(2*Theta[i]) * d_A
            D_alpha = D_alpha/Area
            B_alpha = B_alpha/Area
            Magnitude_of_elongation = np.linalg.norm([[D_alpha, B_alpha], [B_alpha, -D_alpha]]) # Norm of the matrix is the magnitude of elongation
            Elongation_mag.append(Magnitude_of_elongation)
            Mag_limit = 0.5
            
            if view == 'XY':
                pass
            elif view == 'YZ':
                hull = ConvexHull(cell[:,1:3])
                p = cell[hull.vertices,1:3]

            patches = []
            patches.append(Polygon(p, closed=True))
            p = PatchCollection(patches,edgecolor="r",alpha=0.8)
            p.set_color(cpick.to_rgba(Magnitude_of_elongation))
            for k in range(layers):
                if COM_Z < z*(k+1)/layers:
                    ax[k].add_collection(p)
                    break

        if view == 'XY':
            for i in range(layers):
                ax[i].set_xlim([0,size])
                ax[i].set_ylim([0,size])
        elif view == 'YZ':
            for i in range(layers):
                ax[i].set_xlim([0,size])
                ax[i].set_ylim([0,z])
            
        fig.suptitle("Cell Elongation, Time = {}".format(frame_number*1000))
        print('frame_number',frame_number)
        return [ax[i].plot() for i in range(layers)]
    
    def update(frame_number):
        for i in range(layers):
            ax[i].clear()
        frame = th.ReadFrame(inc=nSkip)
        return soft_hard_animation(frame, th, frame_number , filename, size, ax , view , Z) 

    for filename in Files:
        with celldiv.TrajHandle(filename) as th:
            if view == 'XY':
                fig, ax = plt.subplots(1, layers ,  figsize=(5*(layers)+2, 5)) ## only change the xy perspective for now
            elif view == 'YZ':
                fig, ax = plt.subplots(figsize=(size//5,Z//5)) 

            for i in range(layers):
                ax[i].set_xlabel('X')
                ax[i].set_ylabel('Y')
                ax[i].set_aspect('equal', adjustable='box')
            cm1 = mcol.LinearSegmentedColormap.from_list("MyCmapName",['#2A9B2A','r'])
            cnorm = mcol.Normalize(vmin=0,vmax=1)
            cpick = cm.ScalarMappable(norm=cnorm,cmap=cm1)
            cpick.set_array([])
            
            #create an animation
            for _ in range(start_from):
                th.ReadFrame(inc=nSkip)

            ani = animation.FuncAnimation(fig=fig , func=update, frames=range(start_from,th.maxFrames-end,nSkip), interval=50, repeat=False)
            cbar_ax = fig.add_axes([0.93, 0.1, 0.01, 0.8])
            cbar = plt.colorbar(cpick, cax=cbar_ax)
            #cbar = plt.colorbar(cpick, ax=ax[layers-1] , pad =0.1)
            cbar.set_label("$\\Upsilon$", labelpad=10)  # Main label (on the side)
            cbar.ax.set_title('elongated', pad=10, fontsize=12)  # Label above the color bar
            cbar.ax.set_xlabel('Spherical', labelpad=10 , fontsize=12)  # Label below the color bar
            # plt.show()
            # save animation as gif
            if view == 'XY':
                ani.save('elong_colormap_{}.gif'.format(filename[:-4]), writer='Pillow', fps=5)
            elif view == 'YZ':
                ani.save('elong_colormap_YZ_{}.gif'.format(filename[:-4]), writer='Pillow', fps=5)


def Animate_vector_field_elongation_vector(size_z = 28, frame_rate = 5, differenciate_soft_hard = False):
    start_from =0
    end = 0
    size=100
    view = 'XY'

    def visualisation_animation(f,th,frame_number , Filename , size , ax , view = 'XY' , z = 1.9):
        X_COMs = [] # Position vector components - Cell COM
        Y_COMs = []
        U = [] # Elongation vector components - Cell COM
        V = []
        Elongation_mag_1 = []
        Elongation_mag_2 = []

        m = np.vstack(th.cellInd)
        s1 = np.where(m < 0)[0]
        s2 = np.where(m >= 0 )[0]
        # c = [ f[i][:180] for i in range(len(m)) ] #-> Not needed if differeciating by type 
        
        soft = [ f[i][:180] for i in s1 ]
        for cell in soft:
            hull = ConvexHull(cell[:,0:2]) # Finding the boundary points
            Area = hull.volume
            p = cell[hull.vertices,0:2]

            COM = np.mean(p, axis=0)
            X_COMs.append(COM[0])
            Y_COMs.append(COM[1])
            centered_p = p - COM

            # calculations based on paper by Comelles et.al: https://doi.org/10.7554/eLife.57730
            ''' Exx = D_alpha,  Eyy = -D_alpha,  Exy = Eyx = B_alpha '''
            '''D_alpha = 1/A * sum( Cos(2theta) dA )'''
            '''B_alpha = 1/A * sum( Sin(2theta) dA )'''
            D_alpha = np.float64(0)
            B_alpha = np.float64(0)
            Theta = np.arctan2(centered_p[:,1], centered_p[:,0]) # angle of each point with respect to COM
            for i in range(len(Theta)):
               d_A = 0.5 * np.abs(centered_p[i-1][0] * centered_p[i][1] - centered_p[i][0] * centered_p[i-1][1])
               D_alpha += np.cos(2*Theta[i]) * d_A
               B_alpha += np.sin(2*Theta[i]) * d_A
            D_alpha = D_alpha/Area
            B_alpha = B_alpha/Area

            Magnitude_of_elongation = np.linalg.norm([[D_alpha, B_alpha], [B_alpha, -D_alpha]]) 
            Elongation_mag_1.append(Magnitude_of_elongation)
            Mag_of_elongation = np.sqrt(D_alpha**2 + B_alpha**2)
            Elongation_mag_2.append(Mag_of_elongation)
            Elongation_angle = 0.5 * np.arcsin(B_alpha/Magnitude_of_elongation) # angle of elongation vector with respect to x axis
            Elong_angle = np.arctan2(B_alpha, D_alpha)/2
        
            U.append(Mag_of_elongation * np.cos(Elong_angle))
            V.append(Mag_of_elongation * np.sin(Elong_angle))

            # print("Cell index: {}, Elongation magnitude: {}, Elong ang: {}".format(index, Mag_of_elongation, Elong_angle))
            # Plotting the elongation vector field
            patches = []
            patches.append(Polygon(p, closed=True))
            p = PatchCollection(patches, alpha=0.2, edgecolor='blue')
            p.set_color("b") if differenciate_soft_hard else p.set_color(cpick.to_rgba(Magnitude_of_elongation))
            ax.add_collection(p)

        hard = [ f[i][:180] for i in s2 ]
        for cell in hard:
            hull = ConvexHull(cell[:,0:2]) # Finding the boundary points
            Area = hull.volume
            p = cell[hull.vertices,0:2]

            COM = np.mean(p, axis=0)
            X_COMs.append(COM[0])
            Y_COMs.append(COM[1])
            centered_p = p - COM

            # calculations based on paper by Comelles et.al: https://doi.org/10.7554/eLife.57730
            ''' Exx = D_alpha,  Eyy = -D_alpha,  Exy = Eyx = B_alpha '''
            '''D_alpha = 1/A * sum( Cos(2theta) dA )'''
            '''B_alpha = 1/A * sum( Sin(2theta) dA )'''
            D_alpha = np.float64(0)
            B_alpha = np.float64(0)
            Theta = np.arctan2(centered_p[:,1], centered_p[:,0]) # angle of each point with respect to COM
            for i in range(len(Theta)):
               d_A = 0.5 * np.abs(centered_p[i-1][0] * centered_p[i][1] - centered_p[i][0] * centered_p[i-1][1])
               D_alpha += np.cos(2*Theta[i]) * d_A
               B_alpha += np.sin(2*Theta[i]) * d_A
            D_alpha = D_alpha/Area
            B_alpha = B_alpha/Area

            Magnitude_of_elongation = np.linalg.norm([[D_alpha, B_alpha], [B_alpha, -D_alpha]]) 
            Elongation_mag_1.append(Magnitude_of_elongation)
            Mag_of_elongation = np.sqrt(D_alpha**2 + B_alpha**2)
            Elongation_mag_2.append(Mag_of_elongation)
            Elongation_angle = 0.5 * np.arcsin(B_alpha/Magnitude_of_elongation) # angle of elongation vector with respect to x axis
            Elong_angle = np.arctan2(B_alpha, D_alpha)/2
        
            U.append(Mag_of_elongation * np.cos(Elong_angle))
            V.append(Mag_of_elongation * np.sin(Elong_angle))

            # print("Cell index: {}, Elongation magnitude: {}, Elong ang: {}".format(index, Mag_of_elongation, Elong_angle))
            # Plotting the elongation vector field
            patches = []
            patches.append(Polygon(p, closed=True))
            p = PatchCollection(patches,alpha=0.2, edgecolor='red')
            p.set_color("r") if differenciate_soft_hard else p.set_color(cpick.to_rgba(Magnitude_of_elongation))
            ax.add_collection(p)

        print(max(Elongation_mag_1))
        print(max(Elongation_mag_2))
        ax.quiver(X_COMs, Y_COMs, U, V, color = 'k', width = 0.001, pivot = 'mid', headlength=0, headwidth=0, headaxislength=0)
        # plt.plot(X_COMs , Y_COMs)

        if view == 'XY':
            ax.set_xlim([0,size])
            ax.set_ylim([0,size])


        return ax.plot()
    
    def update(frame_number):
        ax.clear()
        frame = th.ReadFrame(inc=nSkip)
        return visualisation_animation(frame, th, frame_number , filename, size, ax , view , size_z)

    for filename in Files:
        with celldiv.TrajHandle(filename) as th:
            if view == 'XY':
                fig, ax = plt.subplots(figsize = (size/5, size/5))
                fig.tight_layout()

            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_aspect('equal', adjustable='box')
            # range of x and y axis
            ax.set_xlim([25,size-25])
            ax.set_ylim([25,size-25])
            # remove the numbers from the axis

            #create an animation
            # for _ in range(start_from):
            #     th.ReadFrame(inc=nSkip)
            ani = animation.FuncAnimation(fig=fig , func=update, frames=range(start_from,th.maxFrames-end,nSkip), interval=50, repeat=False)
            if not differenciate_soft_hard:
                cm1 = mcol.LinearSegmentedColormap.from_list("MyCmapName",['r','b'])
                cnorm = mcol.Normalize(vmin=0,vmax=1)
                cpick = cm.ScalarMappable(norm=cnorm,cmap=cm1)
                cpick.set_array([])
                cbar = plt.colorbar(cpick, ax=ax , pad =0.1 , shrink=0.8)
                cbar.set_label("$\\Upsilon$", labelpad=10)  # Main label (on the side)
                cbar.ax.set_title('Elongated', pad=10, fontsize=12)  # Label above the color bar
                cbar.ax.set_xlabel('Spherical', labelpad=10 , fontsize=12)  # Label below the color bar
            # plt.show()
            #save animation as gif
            if view == 'XY':
                ani.save('elong_colorma_XY_{}.gif'.format(filename[:-4]), writer='Pillow', fps=frame_rate)


def Elongation2(f,frame_number , j , Filename, Plot_distribution = False ,Plot_elongated_cells = False , color_elongated_cells=False , size = 100 , bin_n=50 ): #FIX
    m = np.vstack(th.cellInd)
    c = [ f[i][:180] for i in range(len(m)) ]
    Elongation_mag = []

    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    index_of_elongated_cells = []
    ## color elongation
    cm1 = mcol.LinearSegmentedColormap.from_list("MyCmapName",["r","b"])
    cnorm = mcol.Normalize(vmin=0,vmax=1)
    cpick = cm.ScalarMappable(norm=cnorm,cmap=cm1)
    cpick.set_array([])
    #
    for index , cell in enumerate(c):
        hull = ConvexHull(cell[:,0:2])
        Area = hull.volume
        p = cell[hull.vertices,0:2]
        COM = np.mean(p, axis=0)
        centered_p = p - COM
        '''D_alpha = 1/A * sum( Cos(2theta) dA )'''
        '''B_alpha = 1/A * sum( Sin(2theta) dA )'''
        D_alpha = np.float64(0)
        B_alpha = np.float64(0)
        Theta = np.arctan2(centered_p[:,1], centered_p[:,0])
        for i in range(len(Theta)):
            d_A = 0.5 * np.abs(centered_p[i-1][0] * centered_p[i][1] - centered_p[i][0] * centered_p[i-1][1])
            D_alpha += np.cos(2*Theta[i]) * d_A
            B_alpha += np.sin(2*Theta[i]) * d_A
        D_alpha = D_alpha/Area
        B_alpha = B_alpha/Area

        Magnitude_of_elongation = np.linalg.norm([[D_alpha, B_alpha], [B_alpha, -D_alpha]]) # Norm of the matrix is the magnitude of elongation
        Elongation_mag.append(Magnitude_of_elongation)

        # if index == 3677:
        #     print('Elong for cell 3677',Magnitude_of_elongation)

            
        Mag_limit = 0.5
        
        if Plot_elongated_cells:
            patches = []
            patches.append(Polygon(p, closed=True))
            if Magnitude_of_elongation > Mag_limit:
                p = PatchCollection(patches, color= (1, 128/255, 0 ) ,alpha=0.8)
                ax.add_collection(p)
                index_of_elongated_cells.append(index)
            elif m[index] < 0:
                p = PatchCollection(patches, color= (0.2, 0.6, 1) ,edgecolor="b",alpha=0.8)
                ax.add_collection(p)  
            elif m[index] >=0:
                p = PatchCollection(patches, color= (222/255, 28/255, 0 ) ,alpha=0.8)
                ax.add_collection(p)


        if color_elongated_cells:
            patches = []
            hull = ConvexHull(cell[:,0:2])
            p = cell[hull.vertices,0:2]
            patches.append(Polygon(p, closed=True))
            p = PatchCollection(patches,edgecolor="r",alpha=0.8)
            # for simplex in hull.simplices:
            #     plt.plot(f[mi*192:mi*192 + 180][simplex, 0], f[mi*192:mi*192 + 180][simplex, 1], 'k-')
            p.set_color(cpick.to_rgba(Magnitude_of_elongation))
            ax.add_collection(p)


    # Find the values in array Elongation_mag larger than mag_limit
    Elong_mag_over_limit = [ x for x in Elongation_mag if x > Mag_limit]
    Elongation_mag = np.array(Elongation_mag)
    

    #print(Elong_mag_over_limit)

    if color_elongated_cells:
        ax.set_xlim([0,size])
        ax.set_ylim([0,size])
        plt.colorbar(cpick,label="Area",ax=plt.gca())
        plt.savefig("Elong_colormap_{}_{}.png".format(Filename[:-4], frame_number))

    if Plot_elongated_cells:
        ax.set_xlim([0,size])
        ax.set_ylim([0,size])
        plt.title("Cells with $\\Upsilon > ${}, Time = {}".format(Mag_limit, (frame_number-2)*10000))
        plt.savefig("Elongation_{}_{}.png".format(Filename[:-4], frame_number))
        plt.clf()

    if Plot_distribution:
        plot_distribution_of_elongation(Elongation_mag, bin_n , j , Filename)



Elongation()

#Animate_vector_field_elongation_vector(size_z = 2, frame_rate = 5, differenciate_soft_hard= True)