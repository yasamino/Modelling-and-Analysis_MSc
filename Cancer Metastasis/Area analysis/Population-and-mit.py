import numpy as np

def count_cells(f,th):
    m = np.vstack(th.cellInd)
    total = len(m)
    s1 = np.where(m < 0)[0]
    s2 = np.where(m >= 0 )[0]
    h = [ f[i][:180] for i in s1 ]
    soft = len(h)  
    g = [ f[i][:180] for i in s2 ]
    hard = len(g)
    return soft, hard, total

def Report_mit():
    # Initialize a list to store the total cell counts and corresponding time points
    total_cells = []
    time_points = []
    start=0
    end = 10000
    traj = 10000

    with celldiv.TrajHandle(argv[1]) as th:
        try:
            for frame_number in range(0, th.maxFrames, nSkip):
                if start < frame_number < end:  # Ensure we are within the specified range
                    f = th.ReadFrame(inc=0)  # Read the frame
                    s,h, total = count_cells(f, th)  # Get the total cell count
                    
                    # Store the total count and corresponding time (frame_number)
                    total_cells.append(total)
                    time_points.append(frame_number * traj * 10**(-4))  # Assuming this is the correct time scale
                    
                    # Your existing code continues here...
                    
        except celldiv.IncompleteTrajectoryError:
            print("Stopping...")

    # After the loop, plot the results
    import matplotlib.pyplot as plt

    plt.figure(figsize=(10, 6))
    plt.scatter(time_points, total_cells, label='Total Cells', color='green')

    plt.xlabel('Time')
    plt.ylabel('Cell Count')
    plt.title('Total Cell Count vs Time')
    plt.legend()
    # plt.grid(True)
    plt.tight_layout()

    # Save the plot
    plt.savefig('total_cell_count_vs_time_{}.png'.format(argv[1][:-4]))
    plt.show()
