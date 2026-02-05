def calculate_number_density(f , size, grid , file_name):
    all_coords = np.vstack(f)
    all_coords = np.vstack(all_coords)
    all_com = []
    for mi in range(int(len(all_coords)/180)):
        all_com.append(np.mean(all_coords[mi*180:(mi+1)*180],axis=0))
    all_com = np.array(all_com)[: , 0:2] #2D coordinates
    leng_points = int(size/grid)
    grid_2D = np.zeros(( leng_points , leng_points ))
    for coord in all_com:
        grid_2D[int(coord[0]/grid) , int(coord[1]/grid)]+=1
    Ft = fft2(grid_2D)
    print(Ft)
    IFt = ifft2(Ft/2)
    print(IFt)
    plt.imshow(grid_2D)
    plt.colorbar()
    plt.savefig('grid_{}.png'.format(file_name))
    plt.clf()
    plt.imshow(np.real(IFt))
    plt.colorbar()
    plt.savefig('IFT_{}.png'.format(file_name))
 