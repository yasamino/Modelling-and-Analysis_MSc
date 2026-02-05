import numpy as np
from scipy.spatial import Voronoi

def deMixing_ratio(f,th):
    m = np.vstack(th.cellInd)
    s1 = np.where(m < 0)[0]
    s2 = np.where(m >= 0 )[0]
    h = [ f[i][:180] for i in s1 ]
    len_type1 = len(h)  
    g = [ f[i][:180] for i in s2 ]
    len_type2 = len(g)

    if len_type1 == 0:
        all_coords = np.vstack(g)
    elif len_type2 == 0:
        all_coords = np.vstack(h)
    else:
        all_coords = np.concatenate((h,g), axis=0)

    all_coords = np.vstack(all_coords)
    '''create list of center of mass, form the voronoi diagram, find the neighbors of each cell'''
    all_com = []
    for mi in range(int(len(all_coords)/180)):
        all_com.append(np.mean(all_coords[mi*180:(mi+1)*180],axis=0))
    all_com = np.array(all_com)[: , 0:2]
    vor = Voronoi(all_com)
    '''create list of neighbors for each cell'''
    neighbors = [ [] for _ in range(len(all_com)) ]
    for i in range(len(vor.ridge_points)):
        neighbors[vor.ridge_points[i][0]].append(vor.ridge_points[i][1])
        neighbors[vor.ridge_points[i][1]].append(vor.ridge_points[i][0])

    Dp=0
    for cell in range(len(all_com)):
        Ns = 0
        ''' Check to see if cell index belongs in the range of type 1'''
        if cell < len_type1: 
            #it is a type 1 cell
            for neighbor in neighbors[cell]:
                if neighbor < len_type1:
                    Ns+=1
        else:
            #it is a type 2 cell
            for neighbor in neighbors[cell]:
                if neighbor >= len_type1:
                    Ns+=1

        Nt= len(neighbors[cell]) #total number of neighbors
        Dp+= (2*(Ns/Nt) - 1) #deMixing ratio

    Dp = Dp/len(all_com)
    return Dp



def Animate_correlation_function():
    def Animate_correlation(f, th, frame_number , size):
        m = np.vstack(th.cellInd)
        s1 = np.where(m < 0)[0]
        h = [ f[i][:180] for i in s1 ]
        soft_coords = np.vstack(h)
        soft_com = []
        for soft in range(int(len(soft_coords)/180)):
            soft_com.append(np.mean(soft_coords[soft*180:(soft+1)*180],axis=0))
        soft_com = np.array(soft_com)[: , 0:2] #2D coordinates
        distance = []
        for a in soft_com:
            for b in soft_com:
                distance.append(np.linalg.norm(a-b))
        distance = np.array(distance)
        bins = 100
        r = np.linspace(0, size, bins+1)
        y = [0]*(bins+1)
        for i in distance:
            y[int(i/(size)*bins)]+=1
        N_beta = len(s1)
        rho = N_beta / len(m)

        Area = []
        for cell in h:
            hull = ConvexHull(cell[:,0:2])
            Area.append(np.sqrt(hull.volume))
        print(np.mean(Area))

        plt.title("DP = {} , Time = {}".format( np.round(deMixing_ratio(f,th),4) , (frame_number)*10000))
        plt.xlabel('r')
        plt.ylabel('$g_{ss} (r)$')
        plt.xscale('log')
        return plt.scatter(r[1:],[ y[i]/(2 * np.pi * r[i] * N_beta * rho) for i in range(1,len(y))])


    def update(frame_number):
        ax.clear()
        frame = th.ReadFrame(inc=nSkip)
        return Animate_correlation(frame, th, frame_number , 100)


    for filename in Files:
        with celldiv.TrajHandle(filename) as th:
            fig, ax = plt.subplots()
            #create an animation
            ani = animation.FuncAnimation(fig=fig , func=update, frames=range(0,th.maxFrames,nSkip), interval=30, repeat=False)
            #save animation as gif
            ani.save('Correlation_function_{}.gif'.format(filename[:-4]), writer='Pillow', fps=10)


            
def Animate_Demixing_parameter():
    Dp=[]
    time = []

    def update(frame_number):
        frame = th.ReadFrame(inc=nSkip)
        Dp.append(deMixing_ratio(frame,th))
        time.append(frame_number*10000)
        return plt.plot(time, Dp)
    
    for filename in Files:
        with celldiv.TrajHandle(filename) as th:
            fig, ax = plt.subplots()
            plt.xlabel('Time')
            plt.ylabel('Demixing Parameter')
            plt.ylim(0,1)
            plt.xlim(0,th.maxFrames*10000)
            #create an animation
            ani = animation.FuncAnimation(fig=fig , func=update, frames=range(0,th.maxFrames,nSkip), interval=30, repeat=False)
            #save animation as gif
            ani.save('Demixing_{}.gif'.format(filename[:-4]), writer='Pillow', fps=10)
