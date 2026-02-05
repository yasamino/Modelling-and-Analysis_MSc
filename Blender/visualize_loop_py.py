# A scripts to visualize the files we get with loop_py file
import subprocess
import os



strengths = [ 0.1  ]
epsilon = [ 10 ]
attraction_ranges = [ 0.25]
gamma_friction = [ 0.4 ]
visc_friction = [ 1 ]
repulstion_strengths = [  10 ]
surf_gamma = [ 0.0 ]
#surface_fric = 0.3
lat_pressure = [ 0.4  ]
ang_const = [200 ]
growth_rate = [ 2e-06]
compression = [ 0.21 ]
max_pressure = [ 2.6 ]
angle_constant_relative = [ 1 ]
max_p = 2.6

sigma_array = [1]
rep_range = 0.15
surf = 0.0
#rad = int(input("Radius: "))
rad = 100

for eps in epsilon:
    for rep in repulstion_strengths:
        for visc in visc_friction:
            for gamma in gamma_friction:
                for ang in ang_const:
                    for ang_rel in angle_constant_relative:
                        t=1
                        for gr in growth_rate:
                            for att in strengths:
                                for att_range in attraction_ranges:
                                    for lat in lat_pressure:
                                        for comp in compression:
                                            # if att >1.0 and visc > 0.4:
                                            #     continue
                                            # for surf in surf_gamma:
                                            for rand_scale in [0.05, 0.005]:
                                                # for max_p in max_pressure:
                                                for sigma in sigma_array:
                                                    file_name = "Sep14_Levy_elong_{}_repeat_{}_R{}_2_2.6_{}_1e-3_Jul7_Area_in_pressure_{}_div_stiff0_7_repRange{}_gr{}_ang{}__Att_strength{}_ep{}_sigma{}_attR{}_gamma{}_vis{}_Rep{}_comp{}_lat{}.xyz".format(rand_scale,t, rad, max_p, surf,rep_range, gr, ang, att, eps, sigma, att_range, gamma, visc, rep, comp,lat)
                                                    # -> Copy the file
                                                    # att = att_s / ((att_range) ** 2 - (rep_range)**2 )
                                                    if not (os.path.isfile(file_name)):
                                                        if (os.path.isfile("../bin/" + file_name)):
                                                            command = ['cp', "../bin/" + file_name , './']
                                                            subprocess.run(command)
                                                    # -> Read JSON
                                                    if not (os.path.isdir(file_name[:-4])):
                                                        if os.path.isfile(file_name):
                                                            print("Visualizing file: {}".format(file_name))
                                                            command = ['blender', '--background','CellDiv_Box100_very_large.blend', '--python', 'renderStiff-new.py', '--', file_name] 
                                                            subprocess.run(command)
                                                        else:
                                                            print("File not found: {}".format(file_name))
                                                            continue
                                                    else:
                                                        print("Directory already exists: {}\n".format(file_name[:-4]))
                                                    
                                                    # command = ['python3', 'frame_to_gif.py', file_name]
                                                    # subprocess.run(command)
                                                    #The list of outputs:
                                                    '''
                                                    output_type = input("Enter the type of the output (gif, mp4): ")
                                                    compression = round(float(input("compression (0-1): ")),1 )
                                                    rad = int(input("Radius: "))
                                                    Additional_parameters = input("Additional parameters (optional, press Enter to skip): ")
                                                    if Additional_parameters:
                                                        rep_strength = float(input("Repulsion strength: "))
                                                        att_strength = float(input("Attraction strength: "))
                                                        att_range = float(input("Attraction range: "))
                                                        LJ_sigma = float(input("LJ sigma: "))
                                                        LJ_epsilon = float(input("LJ epsilon: "))
                                                        lat_pressure = float(input("Lateral pressure: "))
                                                        angle_constant = float(input("Angle constant: "))
                                                        growth_rate = float(input("Growth rate: "))
                                                        max_pressure = float(input("Max pressure: "))
                                                        gamma_friction = float(input("Gamma friction: "))
                                                        '''
                                                    

                                                    # Run the subprocess with inputs piped
                                                    if (os.path.isdir(file_name[:-4])):
                                                        if not ( os.path.isfile('animated_' +file_name[:-4] + '.gif')) and not ( os.path.isfile('animated_' + file_name[:-4] + '.mp4')):
                                                            list_of_inputs = ['mp4', comp, rad, True, rep, att, att_range, sigma, eps, lat, ang, gr, max_p, gamma]
                                                            input_str = ''
                                                            for i in range(len(list_of_inputs)):
                                                                input_str += f'{list_of_inputs[i]}\n'
                                                            subprocess.run(['python3', 'frame_to_gif.py', file_name[:-4]], input=input_str, text=True)
                                                        else:
                                                            print("File already exists: animated_" + file_name[:-4] + ".gif or .mp4")
                                                            
                                                


                                                    
                                        
