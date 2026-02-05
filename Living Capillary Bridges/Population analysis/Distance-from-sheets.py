def Animate_distance_from_surface(cutoff_z = 5 , size_z = 30):

    def z_position_report(f,th,frame_number , Filename , size , ax):
        all_coords = np.vstack(f)
        all_coords = np.vstack(all_coords)
        # Just keep the z coordinates
        z_coords = all_coords[:, 2]
        z_com = []
        for mi in range(int(len(z_coords)/180)):
            zCOM=np.mean(z_coords[mi*180:(mi+1)*180],axis=0)
            if (zCOM < cutoff_z) or ((size_z - zCOM )< cutoff_z):
                z_com.append(np.min([zCOM , size_z - zCOM]))
        plt.title("Cell Elongation, Time = {}".format(frame_number*1000))
        plt.xlabel('Distance from the sheets')
        plt.ylabel('Number of cells')
        plt.ylim(0,100)
        plt.xlim(0,cutoff_z)
        #plt.hist(z_com, bins=100)
        return plt.hist(z_com, bins=100)
    
    def update(frame_number):
        ax.clear()
        frame = th.ReadFrame(inc=nSkip)
        return z_position_report(frame, th, frame_number, filename, 100, ax)

    for filename in Files:
        with celldiv.TrajHandle(filename) as th:
            fig, ax = plt.subplots()
            #create an animation
            plt.xlabel('Distance from the sheets')
            plt.ylabel('Number of cells')
            plt.ylim(0,100)
            ani = animation.FuncAnimation(fig=fig , func=update, frames=range(0,th.maxFrames,nSkip), interval=30, repeat=False)
            #save animation as gif
            ani.save('Z_Coords_{}.gif'.format(filename[:-4]), writer='Pillow', fps=10)
