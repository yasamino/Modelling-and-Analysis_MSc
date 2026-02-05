   
def Correlation_function( th , size ,i):
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
    
    plt.scatter(r[1:],[ y[i]/(2 * np.pi * r[i] * N_beta * rho) for i in range(1,len(y))])
    plt.xscale('log')
    plt.savefig('correlation_function_{}.png'.format(i))
            
def Animate_correlation_function(function):
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

        plt.title("DP = {} , Time = {}".format( np.round(function(f,th),4) , (frame_number)*10000))
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
