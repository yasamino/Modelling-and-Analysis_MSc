def Droplet_calculations( view  = "XZ" , filename = Files[0] , Z_precursor = 1.5 , Z_contact = 10.5 , Size_Box = 100 ):
    ''' Input 
    view: "XZ" or "YZ" - the view of the droplet from the sides
    Z_ precursor: The height at which the precursor will be located, close to the surface
    Z_contact: The height at which the contact line will be located, close to the surface and higher than the prescursor
    Size_Box: Size of box 
    '''
    with celldiv.TrajHandle(filename) as th:
        frameCount = 0
        mean_precursor_radius = []
        mean_contact_radius = []
        mean_max_precursor_radii = []
        mean_max_contact_radii = []

        unit_conversion = 7/(0.9) # 7 micro meters = 0.9 sim unit
        


        try:
            for i in range(int(th.maxFrames/nSkip)+1): # i for each frame written to file
                frameCount += 1
                if frameCount > args.num_frames:
                    break
                f=th.ReadFrame(inc=nSkip ) # f is the list of all the cells in the frame
                # f_count+= len(f)

                m = np.vstack(th.cellInd)
                all_coords = np.vstack(f)
                all_coords = np.vstack(all_coords)
                all_com = []
                for mi in range(int(len(all_coords)/180)):
                    all_com.append(np.mean(all_coords[mi*180:(mi+1)*180],axis=0))
                all_com = np.array(all_com)# * unit_conversion

                precursor_slice = [ all_com[i]-(Size_Box)/2 for i in range(len(all_com))  if all_com[i][2] < Z_precursor ] # get the cells that are below the precursor height
                mean_precursor_radius.append(np.mean(np.sqrt([np.linalg.norm(precursor_slice[i][0:2])**2 for i in range(len(precursor_slice))])) )

                # precursor_radius_max = np.max([np.linalg.norm(precursor_slice[i][0:2]) for i in range(len(precursor_slice))])
                contact_slice = [ all_com[i]-(Size_Box)/2 for i in range(len(all_com)) if (Z_contact-1 < all_com[i][2] and all_com[i][2] < Z_contact) ] # get the cells that are below the contact height
                mean_contact_radius.append(np.mean(np.sqrt([np.linalg.norm(contact_slice[i][0:2])**2 for i in range(len(contact_slice))])))
                # Contact_radius_max = np.max([np.linalg.norm(contact_slice[i][0:2]) for i in range(len(contact_slice))])
                # print("Contact radius: ", Contact_radius)
                
                mean_max_contact_radii.append(np.mean(np.sort(np.sqrt([np.linalg.norm(contact_slice[i][0:2])**2 for i in range(len(contact_slice))]))[-20:-10]))
                mean_max_precursor_radii.append(np.mean(np.sort(np.sqrt([np.linalg.norm(precursor_slice[i][0:2])**2 for i in range(len(precursor_slice))]))[-20:-10]))



                if len(args.inds) > 0:
                    pass

        except celldiv.IncompleteTrajectoryError:
            print ("Stopping...") 

    fig, ax = plt.subplots(figsize=(10, 6))
    # convert to simulation units
    mean_precursor_radius = [i*unit_conversion for i in mean_precursor_radius]
    # max_precursor_radius = [i*unit_conversion for i in max_precursor_radius]
    mean_contact_radius = [i*unit_conversion for i in mean_contact_radius]
    # max_contact_radius = [i*unit_conversion for i in max_contact_radius]
    mean_max_precursor_radii = [i*unit_conversion for i in mean_max_precursor_radii]
    mean_max_contact_radii = [i*unit_conversion for i in mean_max_contact_radii]

    # find the average of the 10 largest values in contact radius
    



    color_precursor = 'k'#"#1f77b4"  # Blue for precursor
    color_contact = 'k' #"#ff7f0e"    # Orange for contact

    # X-axis: frame numbers
    frames = list(range(len(mean_precursor_radius)))

    # Plot: Linear scale
    fig, ax = plt.subplots()
    ax.plot(frames[1:], mean_precursor_radius[1:], label='Mean Precursor Radius', marker='o', color=color_precursor)
    ax.plot(frames[1:], mean_max_precursor_radii[1:] , label='Max Precursor Radius(avg10)', marker='o', linestyle='--', color=color_precursor)
    ax.plot(frames[1:], mean_contact_radius[1:], label='Mean Contact Radius', marker='<', color=color_contact)
    ax.plot(frames[1:], mean_max_contact_radii[1:] , label='Max Contact Radius(avg10)', marker='<', linestyle='--', color=color_contact)
    ax.set_xlabel('Time (frames)')
    ax.set_ylabel('Radius [$\\mu$m]')
    ax.set_title('Droplet Radius Calculations (Linear Scale)')
    ax.legend()
    plt.tight_layout()
    plt.savefig('Droplet_Radius_Calculations_{}_linear.svg'.format(filename[:-4]))

    # Plot: Log-log scale
    fig, ax = plt.subplots()
    ax.plot(frames[1:], mean_precursor_radius[1:], label='Mean Precursor Radius', marker='o', color=color_precursor)
    ax.plot(frames[1:], mean_max_precursor_radii[1:], label='Max Precursor Radius', marker='o', linestyle='--', color=color_precursor)
    ax.plot(frames[1:], mean_contact_radius[1:], label='Mean Contact Radius', marker='<', color=color_contact)
    ax.plot(frames[1:], mean_max_contact_radii[1:], label='Max Contact Radius', marker='<', linestyle='--', color=color_contact)
    ax.set_xlabel('Time (frames)')
    ax.set_ylabel('Radius [$\\mu$m]')
    ax.set_title('Droplet Radius Calculations (Log-Log Scale)')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend()
    plt.tight_layout()
    plt.savefig('Droplet_Radius_Calculations_{}_loglog.png'.format(filename[:-4]))

    # Save data to CSV
    data = {
        'Frame Number': frames,
        'Mean Precursor Radius': mean_precursor_radius,
        'Max Precursor Radius': mean_max_precursor_radii,
        'Mean Contact Radius': mean_contact_radius,
        'Max Contact Radius': mean_max_contact_radii
    }
    df = pd.DataFrame(data)
    df.to_csv('Droplet_Radius_Calculations_{}.csv'.format(filename[:-4]), index=False)





    # we first start by reading the trajectory file and extracting the frames
    # 
