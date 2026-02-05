

def Report_count():
    # Initialize a list to store the total cell counts and corresponding time points
    total_cells = []
    total_new=[]
    total_old=[]
    time_points = []
    start=2
    end = 100000000
    traj = 10000
    count_interval = 1
    initial_count =0

    with celldiv.TrajHandle(argv[1]) as th:
        try:
            for frame_number in range(0, th.maxFrames, nSkip):
                # frame = th.ReadFrame(inc=nSkip)
                print(th.maxFrames)
                # if start <= frame_number <= end:  # Ensure we are within the specified range
                # if frame_number % count_interval == 0: 
                print(frame_number)
                f = th.ReadFrame(inc=1)  # Read the frame
                if frame_number == start:
                    initial_count = len(f)
                
                # s,h, total = count_cells(f, th)  # Get the total cell count
                all,total, total_n, total_o = count_surface_cells(f, th, initial_count)  # Get the total cell count on the surface
                
                # Store the total count and corresponding time (frame_number)
                total_cells.append(total)
                total_new.append(total_n)
                total_old.append(total_o)
                time_points.append(frame_number)  # Assuming this is the correct time scale
                    
                        # Your existing code continues here...
                    
        except celldiv.IncompleteTrajectoryError:
            print("Stopping...")

    # After the loop, plot the results
    import matplotlib.pyplot as plt

    plt.figure(figsize=(10, 6))
    plt.scatter(time_points, total_cells, label='Total Cells', color='black')
    plt.scatter(time_points, total_new, label='New Cells', color='orange')
    plt.scatter(time_points, total_old, label='Old Cells', color='green')
    #fit line to the data
    # z = np.polyfit(time_points[3:], total_cells[3:], 1) 
    # p = np.poly1d(z)
    # plt.plot(time_points[3:],p(time_points[3:]),"k--", label='Fit: y={:.2f}x+{:.2f}'.format(z[0],z[1]))

    plt.xlabel('Time')
    plt.ylabel('Cell Count')
    plt.title('surface Cell Count vs Time')
    plt.legend()
    # plt.grid(True)
    plt.tight_layout()

    # Save the plot
    plt.savefig('surface_cell_count_vs_time_{}.png'.format(argv[1][:-4]))
    plt.show()



def Report_count_Bulk():
    # Initialize a list to store the total cell counts and corresponding time points
    total_bulk=[]

    time_points = []
    start=2
    end = 100000000
    traj = 10000
    count_interval = 1
    initial_count =0

    with celldiv.TrajHandle(argv[1]) as th:
        try:
            for frame_number in range(0, th.maxFrames, nSkip):
                # frame = th.ReadFrame(inc=nSkip)
                print(th.maxFrames)
                # if start <= frame_number <= end:  # Ensure we are within the specified range
                # if frame_number % count_interval == 0: 
                print(frame_number)
                f = th.ReadFrame(inc=1)  # Read the frame
                if frame_number == start:
                    initial_count = len(f)
                
                # s,h, total = count_cells(f, th)  # Get the total cell count
                all,total, total_n, total_o = count_surface_cells(f, th, initial_count)  # Get the total cell count on the surface
                
                # Store the total count and corresponding time (frame_number)
                total_bulk.append(all-total)
                
                time_points.append(frame_number)  # Assuming this is the correct time scale
                    
                        # Your existing code continues here...
                    
        except celldiv.IncompleteTrajectoryError:
            print("Stopping...")

    # After the loop, plot the results
    import matplotlib.pyplot as plt

    plt.figure(figsize=(10, 6))
    plt.scatter(time_points, total_bulk, label='Total Cells in Bulk', color='black')
    #fit line to the data
    # z = np.polyfit(time_points[3:], total_cells[3:], 1) 
    # p = np.poly1d(z)
    # plt.plot(time_points[3:],p(time_points[3:]),"k--", label='Fit: y={:.2f}x+{:.2f}'.format(z[0],z[1]))

    plt.xlabel('Time')
    plt.ylabel('Cell Count')
    plt.title('Cell Count vs Time')
    plt.legend()
    # plt.grid(True)
    plt.tight_layout()

    # Save the plot
    plt.savefig('Bulk_cell_count_vs_time_{}.png'.format(argv[1][:-4]))
    plt.show()

 