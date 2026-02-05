def count_surface_cells(f,th, initial_count):
    # f = th.ReadFrame(inc=nSkip )
    m = np.vstack(th.cellInd)
    c = [ f[i][:180] for i in range(len(m)) ]
    heights = [np.mean(cell, axis=0)[2] for cell in c]
    
    z_min = min(heights) + 1
    z_max = max(heights) - 1
    total_surface = sum(1 for cell in c if (np.mean(cell, axis=0)[2] > z_max or 
                                    np.mean(cell, axis=0)[2] < z_min))
    total_surface_old = sum(1 for cell in c[:initial_count] if (np.mean(cell, axis=0)[2] > z_max or 
                                    np.mean(cell, axis=0)[2] < z_min))
    total_surface_new = total_surface - total_surface_old
    print("Min and Max height of cells: ", np.min(heights), np.max(heights))
    total = len(m)

    return total , total_surface, total_surface_new, total_surface_old

