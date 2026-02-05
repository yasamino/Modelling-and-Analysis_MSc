import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import ConvexHull
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import alphashape # type: ignore
from builtins import len



def Find_boundary_points(f,f_count,size , plot_or_not , m):
    fig = plt.figure(1 , figsize=(5,5))
    ax = fig.add_subplot()
    s1 = np.where(m >= 0)[0]  #( m>0 for hard cells , m<0 for soft cells)
    h = [ f[i][:180] for i in s1 ] #coordinates of nodes on each soft cell
    soft_outer_boundary = []
    count=0
    for cell in h:
        patches=[]
        count+=1
        hull = ConvexHull(cell[:,0:2])
        p = cell[hull.vertices,0:2]
        soft_outer_boundary.extend(p)
        patches.append(Polygon(p, closed=True))
        p = PatchCollection(patches, color = "b" ,edgecolor="b",alpha=0.8)
        ax.add_collection(p)

    ''' Perform concave hull'''
    alpha_shape = alphashape.alphashape(soft_outer_boundary, 1.2) 
    arr = alpha_shape.exterior.coords.xy
    if plot_or_not:
        plt.xlim(0,size)
        plt.ylim(0,size)
        plt.plot(arr[0], arr[1], 'r-')
        plt.savefig('boundary_points{}_clf.svg'.format(f_count))
        plt.clf()
    
    return arr


def Calculate_line_Roughness(f,f_count,size , l ):
    arr = Find_boundary_points(f,f_count,size , plot_or_not = False)
    L0 = 2* size # Total length of line , two sides

    x = arr[0]
    y = arr[1]
    
    segment_number =0
    segment_array = []
    while segment_number < len(x):
        segment_array.append([ [x[segment_number-1],y[segment_number-1]] ])
        L=0
        while L < l :
            L += np.abs(y[segment_number] - y[segment_number-1]) # the vertical distance between two lines
            segment_array[-1].append([x[segment_number],y[segment_number]])
            segment_number+=1
    # Now we have the segments of the line
    # Calculate the roughness of each segment
    Roughness = []
    for segment in segment_array:
        x_seg = np.array([i[0] for i in segment])
        y_seg = np.array([i[1] for i in segment])
        mean_h = np.mean(np.abs(x_seg - size/2)) # mean of distance of x to center
        rough_2 = np.square(x_seg  - mean_h)
        Roughness.append(np.mean(rough_2))
    return np.mean(np.sqrt(Roughness))

def Calculate_line_Roughness_var(f,f_count,size , l ):
    arr = Find_boundary_points(f,f_count,size , plot_or_not = False)
    L0 = 2* size # Total length of line , two sides

    x = arr[0]
    y = arr[1]
    
    segment_number =0
    segment_array = []
    while segment_number < len(x):
        segment_array.append([ [x[segment_number-1],y[segment_number-1]] ])
        L=0
        while L < l :
            L += np.abs(y[segment_number] - y[segment_number-1]) # the vertical distance between two lines
            segment_array[-1].append([x[segment_number],y[segment_number]])
            segment_number+=1
    # Now we have the segments of the line
    # Calculate the roughness of each segment
    Roughness = []
    for segment in segment_array:
        x_seg = np.array([i[0] for i in segment])
        y_seg = np.array([i[1] for i in segment])
        mean_h = np.mean(np.abs(x_seg - size/2)) # mean of distance of x to center
        rough_2 = np.square(x_seg  - mean_h)
        Roughness.append(np.mean(rough_2))
    return np.mean(np.sqrt(Roughness))

'''Scatter plot the boundary of the cells and save the figures so they can be analysed using 
   canny edge detection in File SurfaceAnalysis.ipynb'''

def Plt_boundary_mahmood(f,f_count ,m):
    fig = plt.figure(1 , figsize=(5,5))
    ax = fig.add_subplot()
    s1 = np.where(m < 0)[0]
    h = [ f[i][:180] for i in s1 ] #coordinates of nodes on each soft cell
    soft_outer_boundary = []
    soft_outer_boundary = np.vstack(h)
    plt.plot(soft_outer_boundary[: , 0] , soft_outer_boundary[: , 1] , '.')
    plt.xlim(0,100)
    plt.ylim(0,100)
    plt.axis('off')
    plt.savefig('boundary_mah_{}.png'.format(f_count))
    plt.clf()
