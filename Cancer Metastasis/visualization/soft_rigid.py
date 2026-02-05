
def Animate_soft_hard_layers(): #interactive - change the parameters
    start_from = 0
    end = 0
    
    Z = 50#float(input("Size of the box - Z axis : \n"))
    layers = int(input("Number of layers : \n"))
    size = 100#int(input("Size of the box - X and Y : \n"))
    view = 'XY'#input("view : \n")
    ratio = 0.25#float(input("Ratio of output: \n"))

    def soft_hard_animation(f,th,frame_number , Filename , size , ax , view = 'XY' , z = 1.9):
        m = np.vstack(th.cellInd)
        s1 = np.where(m < 0)[0]
        s2 = np.where(m >= 0 )[0]

        h = [ f[i][:180] for i in s1 ]
        for cell in h:
            COM_Z = np.mean(cell[:,2], axis=0) 
            patches=[]
            hull = ConvexHull(cell[:,0:2])
            p = cell[hull.vertices,0:2]
            patches.append(Polygon(p, closed=True))
            p = PatchCollection(patches, color = (0.2, 0.6, 1) ,edgecolor="b",alpha=0.8)
            ax[int(COM_Z/Z * layers)].add_collection(p)

        g = [ f[i][:180] for i in s2 ]
        for cell in g:
            COM_Z = np.mean(cell[:,2], axis=0) 
            patches=[]
            hull = ConvexHull(cell[:,0:2])
            p = cell[hull.vertices,0:2]
            patches.append(Polygon(p, closed=True))
            p = PatchCollection(patches, color= (0.2, 0.6, 1) ,edgecolor="b",alpha=0.8)##(222/255, 28/255, 0 ) ,alpha=0.8)
            ax[int(COM_Z/Z * layers)].add_collection(p)

        if view == 'XY':
            for i in range(layers):
                ax[i].set_xlim([size*(1-ratio)/2,size*(1+ratio)/2])
                ax[i].set_ylim([size*(1-ratio)/2,size*(1+ratio)/2])

        elif view == 'YZ':
            for i in range(layers):
                ax[i].set_xlim([0,size])
                ax[i].set_ylim([0,z])
            
        fig.suptitle("Boundary cells, Time = {}".format(frame_number*1000))
        print('frame_number',frame_number)
        return [ax[i].plot() for i in range(layers)]
    
    def update(frame_number):
        for i in range(layers):
            ax[i].clear()
            ax[i].set_xlabel('X')
            ax[i].set_ylabel('Y')
            ax[i].set_title('Layer = {}'.format(i+1))
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

            #create an animation
            for _ in range(start_from):
                th.ReadFrame(inc=nSkip)

            ani = animation.FuncAnimation(fig=fig , func=update, frames=range(start_from,th.maxFrames-end,nSkip), interval=50, repeat=False)
            if view == 'XY':
                ani.save('boundary_{}.gif'.format(filename[:-4]), writer='Pillow', fps=5)
            elif view == 'YZ':
                ani.save('boundary_YZ_{}.gif'.format(filename[:-4]), writer='Pillow', fps=5)


def color_soft_hard(f,th,frame_number , Filename , size):
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    Area_soft=[]
    Area_hard=[]

    '''cell node coordinates soft (h) - hard(g)'''
    m = np.vstack(th.cellInd)
    s1 = np.where(m < 0)[0]
    s2 = np.where(m >= 0 )[0]

    h = [ f[i][:180] for i in s1 ]
    for cell in h:
        patches=[]
        hull = ConvexHull(cell[:,0:2])
        p = cell[hull.vertices,0:2]
        patches.append(Polygon(p, closed=True))
        p = PatchCollection(patches, color = (0.2, 0.6, 1) ,edgecolor="b",alpha=0.8)
        ax.add_collection(p)
        Area_soft.append(hull.volume)
    g = [ f[i][:180] for i in s2 ]
    for cell in g:
        patches=[]
        hull = ConvexHull(cell[:,0:2])
        p = cell[hull.vertices,0:2]
        patches.append(Polygon(p, closed=True))
        p = PatchCollection(patches, color= (222/255, 28/255, 0 ) ,alpha=0.8)
        ax.add_collection(p)
        Area_hard.append(hull.volume)

    ax.set_xlim([0,size])
    ax.set_ylim([0,size])
    plt.title("DC = {}% , Time = {}".format(int( (sum(Area_soft) + sum(Area_hard) )/(size*size)*100) , (frame_number-2)*10000))
    plt.savefig("{}_{}.png".format(Filename, frame_number))
    plt.clf()

 