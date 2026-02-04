import matplotlib as plt
from scipy.spatial import ConvexHull
import matplotlib.colors as mcol
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon
from matplotlib import cm
import numpy as np

def Color_Area_dist(f,th, size): 
    fig = plt.figure()
    ax = fig.add_subplot()
    Area = []
    cm1 = mcol.LinearSegmentedColormap.from_list("MyCmapName",["r","b"])
    cnorm = mcol.Normalize(vmin=0,vmax=3)
    cpick = cm.ScalarMappable(norm=cnorm,cmap=cm1)
    cpick.set_array([])
    m = np.vstack(th.cellInd)
    c = [ f[i][:180] for i in range(len(m)) ]
    for cell in c:
        patches=[]
        hull = ConvexHull(cell[:,0:2])
        p = cell[hull.vertices,0:2]
        patches.append(Polygon(p, closed=True))
        p = PatchCollection(patches,edgecolor="r",alpha=0.8)
        # for simplex in hull.simplices:
        #     plt.plot(f[mi*192:mi*192 + 180][simplex, 0], f[mi*192:mi*192 + 180][simplex, 1], 'k-')
        Area.append(hull.volume) #returns the area for 2D case
        p.set_color(cpick.to_rgba(hull.volume))
        ax.add_collection(p)
    ax.set_xlim([0,size])
    ax.set_ylim([0,size])
    plt.colorbar(cpick,label="Area",ax=plt.gca())
    plt.title("Degree of confluence $\\rho$ = {}%".format(int(sum(Area)/(40*40)*100)))
    fig.savefig('Area distribution for rho = {}.svg'.format(int(sum(Area)/(40*40)*100)))