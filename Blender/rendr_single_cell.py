import bpy
import csv
import sys
import os
import argparse
import numpy as np


import sys
import subprocess

user_site = os.path.expanduser("~/.local/lib/python3.11/site-packages")
if user_site not in sys.path:
    sys.path.append(user_site)

# Blender Snap-safe site-packages path
bl_ver = "4.5.1"
py_ver = "3.11"  # match the folder name inside python/lib
site_path = os.path.expanduser(f"~/.config/blender/{bl_ver}/python/lib/python{py_ver}/site-packages")

if site_path not in sys.path:
    sys.path.append(site_path)


try:
    import pandas as pd
except ImportError:
    python_exe = sys.executable  # Blender's own Python interpreter
    subprocess.check_call([python_exe, "-m", "ensurepip", "--upgrade"])
    subprocess.check_call([python_exe, "-m", "pip", "install", "--upgrade", "pip"])
    subprocess.check_call([python_exe, "-m", "pip", "install", "pandas"])
    bl_ver = "4.5.1"
    py_ver = "3.11"  # match the folder name inside python/lib
    site_path = os.path.expanduser(f"~/.config/blender/{bl_ver}/python/lib/python{py_ver}/site-packages")

    if site_path not in sys.path:
        sys.path.append(site_path)

    import pandas as pd


import pandas as pd


# Try to import matplotlib; if it fails, install it into Blender's Python
try:
    import matplotlib
except ImportError:
    python_exe = sys.executable  # Blender's own Python interpreter
    subprocess.check_call([python_exe, "-m", "ensurepip", "--upgrade"])
    subprocess.check_call([python_exe, "-m", "pip", "install", "--upgrade", "pip"])
    subprocess.check_call([python_exe, "-m", "pip", "install", "matplotlib"])
    bl_ver = "4.5.1"
    py_ver = "3.11"  # match the folder name inside python/lib
    site_path = os.path.expanduser(f"~/.config/blender/{bl_ver}/python/lib/python{py_ver}/site-packages")

    if site_path not in sys.path:
        sys.path.append(site_path)

    import matplotlib as mpl


import matplotlib as mpl
import matplotlib.pyplot as plt
sys.path.append("/home/yasamin/multigpu/CellSim3D/scripts")
import celldiv


argv = sys.argv

if "--" not in argv:
    print("ERROR: No arguments provided to script")
    sys.exit(80)
else:
    a = argv.index("--")
    argv = argv[a + 1:]


helpString = """
Run as:
blender --background --python %s --
[options]
""" % __file__

parser = argparse.ArgumentParser(description=helpString)

parser.add_argument("trajPath", type=str,
                    help="Trajectory path. Absolute or relative.")

parser.add_argument("-s", "--smooth", action='store_true',
                    help="Do smoothing (really expensive and doesn't look as good)")

parser.add_argument("-k", "--skip", type=int, required=False,
                    help="Trajectory frame skip rate. E.g. SKIP=10 will only \
                    render every 10th frame.",
                    default=1)

parser.add_argument("-nc", "--noclear", type=bool, required=False,
                    help="specifying this will not clear the destination directory\
                    and restart rendering.",
                    default=False)

parser.add_argument("--min-cells", type=int, required=False,
                    help='Start rendering when system has at least this many cells',
                    default=1)

parser.add_argument("--inds", type=int, required=False, nargs='+',
                    help="Only render cells with these indices",
                    default=[])

parser.add_argument("-nf", "--num-frames", type=int, required=False,
                    help="Only render this many frames.",
                    default=sys.maxsize)

parser.add_argument("-r", "--res", type=int, default=1, required=False,
                    help='Renders images with resolution RES*1080p. RES>=1. \
                    Use 2 for 4k. A high number will devour your RAM.')

parser.add_argument("-cc", "--cell-color", type=int, nargs=3, required=False,
                    default=[82, 38, 123],
                    help="RGB values of cell color. From 0 to 255")

parser.add_argument("-bc", "--background-color", type=int, nargs=3,
                    required=False, default=[255,255,255],
                    help="RGB values of cell color. From 0 to 255")

parser.add_argument("-si", "--specular-intensity", type=float, required=False,
                    default = 0.0,
                    help="Set specular-intensity (shininess). From 0.0 to 1.0")

args = parser.parse_args(argv)

argv = sys.argv
imageindex = 0
firstfaces = []

# Get the active world
world = bpy.context.scene.world

# Set the background color
world.use_nodes = True
# bg_node = world.node_tree.nodes['Background']
# bg_node.inputs['Color'].default_value = [(1.0/255.0)*c for c in args.background_color] + [1.0]


# Make render transparent
bpy.context.scene.render.film_transparent = True
bpy.context.scene.render.image_settings.file_format = 'PNG'
bpy.context.scene.render.image_settings.color_mode = 'RGBA'

bpy.context.view_layer.use_sky = True



doSmooth = args.smooth
if doSmooth:
    print("Doing smoothing. Consider avoiding this feature...")


if (args.res < 1):
    print("ERROR: invalid resolution factor")
    sys.exit()

bpy.data.scenes["Scene"].render.resolution_x*=args.res
bpy.data.scenes["Scene"].render.resolution_y*=args.res

#scene = bpy.context.scene
#scene.render.film_transparent = True  # Enable transparent background
#scene.render.image_settings.file_format = 'PNG'
#scene.render.image_settings.color_mode = 'RGBA'  # Must include Alpha

with open('C180_pentahexa.csv', newline='') as g:
    readerfaces = csv.reader(g, delimiter=',')
    for row in readerfaces:
        firstfaces.append([int(v) for v in row])
        

filename = os.path.realpath(args.trajPath)
basename = "./singlecells"+ os.path.splitext(filename)[0] + "/images/CellDiv_"

nSkip = args.skip

if nSkip > 1:
    print("Rendering every %dth" % nSkip, "frame...")


noClear = args.noclear

sPath = os.path.splitext(filename)[0] + "/images/"

if not noClear and os.path.exists(sPath):
    for f in os.listdir(sPath):
        os.remove(sPath+f)

cellInds = []
minInd = args.min_cells - 1
if len(args.inds) > 0:
    minInd = max(minInd, min(args.inds))

stopAt = args.num_frames

# Set material color
bpy.data.materials['Material'].diffuse_color = [ (1/255.0) * c for c in args.cell_color] + [1.0]
bpy.data.materials['Material'].specular_intensity = args.specular_intensity



def locate_and_plot_by_shape_index(df_F, filename, bins, samples=10, nSkip=1):
    f_count = 0
    base_filename = os.path.splitext(os.path.basename(filename))[0]
    cmap = plt.get_cmap('coolwarm')
    norm = mpl.colors.Normalize(vmin=4.5, vmax=10.5)

    with celldiv.TrajHandle(filename) as th:
        frameCount = 1
        try:
            for i in range(int(th.maxFrames/nSkip)+1):  # i for each frame written to file
                f = th.ReadFrame(inc=nSkip)
                print('frame:', frameCount)

                m = np.vstack(th.cellInd)
                f = [f[j][:180] for j in range(len(m))]
                n_cells = len(f)

                X_com = np.mean([np.mean(cell[:, 0]) for cell in f])
                Y_com = np.mean([np.mean(cell[:, 1]) for cell in f])

                area_vals = df_F["Area"][f_count:f_count + n_cells]
                volume_vals = df_F["Volume"][f_count:f_count + n_cells]

                with np.errstate(divide='ignore', invalid='ignore'):
                    shape_index = np.where(volume_vals > 0, area_vals / (volume_vals ** (2 / 3)), np.nan)
                shape_index = np.where(np.isfinite(shape_index), shape_index, np.nan)

                count = np.zeros(7) # first 4-5 , second 5-6, third 6-7, fourth 7-8, fifth 8-9, sixth 9-10

                for index, cell in enumerate(f):
                    si = shape_index[index]
                    if not np.isfinite(si):
                        continue

                    com = np.mean(cell, axis=0)
                    cell[:, 0] -= com[0]
                    cell[:, 1] -= com[1] 
                    cell[:, 2] -= com[2] 

                    if count[ int(si-4) ] < 2:
                        fig = plt.figure()
                        ax = fig.add_subplot(111, projection='3d')
                        ax.scatter(cell[:, 0], cell[:, 1], cell[:, 2], color =cmap(norm(si)), cmap=cmap, norm=norm, s=10)
                        # print('nodes: ',len(cell[:,0]))
                        ax.set_box_aspect([1, 1, 1])  # Equal aspect ratio
                        ax.axis('off')
                        # print(area_vals, volume_vals)
                        print(f"Cell {index} has shape index {si:.2f}, area {area_vals.iloc[index]:.2f}, volume {volume_vals.iloc[index]:.2f}")
                        plt.title(f"Shape Index {si:.2f}, Area {area_vals.iloc[index]:.2f}, Volume {volume_vals.iloc[index]:.2f}")
                        filename_out = f"{np.round(si,2)}_Frame{frameCount}_{base_filename}.svg"
                        # Add colorbar
                        sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
                        sm.set_array([])  # required for matplotlib < 3.1
                        cbar = plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
                        cbar.set_label("Shape Index")
                        plt.savefig(filename_out)
                        plt.close(fig)
                        plotting_in_blender(cell,si,index, cmap(norm(si)),i )
                        count[ int(si-4) ] += 1


                    if all(co == 2 for co in count):
                        break

                f_count += len(f)
                frameCount += 1

        except celldiv.IncompleteTrajectoryError:
            print("Trajectory incomplete. Ending early.")


# def plotting_in_blender(verts_np, si, mi, color_rgba, frameCount, dot_radius=0.02):
#     # --- Convert vertices ---
#     verts = [tuple(map(float, v)) for v in np.asarray(verts_np).reshape(-1, 3)]
#     faces1 = [tuple(int(v) for v in row) for row in firstfaces]

#     # If faces are 1-based, shift to 0-based
#     if max(max(f) for f in faces1) == len(verts):
#         faces1 = [tuple(v-1 for v in f) for f in faces1]

#     # Sanity check
#     if max(max(f) for f in faces1) >= len(verts):
#         raise ValueError("Face indices out of range for given vertices.")

#     # --- Material (shared for both surface and dots) ---
#     rgba = list(color_rgba) if len(color_rgba) == 4 else list(color_rgba) + [1.0]
#     mat = bpy.data.materials.new(name=f"CellMat_{mi}_{frameCount}")
#     mat.use_nodes = True
#     bsdf = mat.node_tree.nodes.get("Principled BSDF") or mat.node_tree.nodes.new("ShaderNodeBsdfPrincipled")
#     bsdf.inputs["Base Color"].default_value = rgba
#     bsdf.inputs["Alpha"].default_value = 0.4
#     mat.blend_method  = 'BLEND'
#     mat.shadow_method = 'HASHED'

#     # --- Surface mesh ---
#     surf_mesh = bpy.data.meshes.new(f'cellMesh_{mi}_{frameCount}')
#     surf_mesh.from_pydata(verts, [], faces1)
#     surf_mesh.validate(clean_customdata=True)
#     surf_mesh.update()

#     surf_obj = bpy.data.objects.new(f'cellObject_{mi}_{frameCount}', surf_mesh)
#     surf_obj.data.materials.append(mat)
#     bpy.context.scene.collection.objects.link(surf_obj)

#     # --- Points mesh (no faces) ---
#     pts_mesh = bpy.data.meshes.new(f'pointsMesh_{mi}_{frameCount}')
#     pts_mesh.from_pydata(verts, [], [])
#     pts_mesh.update()

#     pts_obj = bpy.data.objects.new(f'pointsObject_{mi}_{frameCount}', pts_mesh)
#     bpy.context.scene.collection.objects.link(pts_obj)

#     # --- Sphere to instance on points ---
#     bpy.ops.mesh.primitive_ico_sphere_add(radius=dot_radius, location=(0, 0, 0))
#     ico = bpy.context.active_object
#     ico.name = f"dot_proto_{mi}_{frameCount}"
#     ico.data.materials.append(mat)

#     # --- Geometry Nodes: instance sphere on each vertex ---
#     gn = bpy.data.node_groups.new(f"GN_InstDots_{mi}_{frameCount}", 'GeometryNodeTree')
#     mod = pts_obj.modifiers.new(name="GN_InstDots", type='NODES')
#     mod.node_group = gn

#     nodes = gn.nodes
#     links = gn.links
#     nodes.clear()

#     ng_input  = nodes.new('NodeGroupInput')
#     ng_output = nodes.new('NodeGroupOutput')
#     gn.inputs.new('NodeSocketGeometry', 'Geometry')
#     gn.outputs.new('NodeSocketGeometry', 'Geometry')

#     obj_info = nodes.new('GeometryNodeObjectInfo')
#     obj_info.inputs['As Instance'].default_value = True
#     obj_info.inputs['Object'].default_value = ico

#     inst_on_points = nodes.new('GeometryNodeInstanceOnPoints')

#     links.new(ng_input.outputs['Geometry'], inst_on_points.inputs['Points'])
#     links.new(obj_info.outputs['Geometry'], inst_on_points.inputs['Instance'])
#     links.new(inst_on_points.outputs['Instances'], ng_output.inputs['Geometry'])

#     # --- Render ---
#     imagename = basename + f"Figure_number_{frameCount}_{si:.2f}.png"
#     bpy.context.scene.render.filepath = imagename
#     bpy.ops.render.render(write_still=True)

#     # --- Cleanup ---
#     bpy.data.objects.remove(surf_obj, do_unlink=True)
#     bpy.data.objects.remove(pts_obj, do_unlink=True)
#     bpy.data.objects.remove(ico, do_unlink=True)
#     bpy.data.meshes.remove(surf_mesh, do_unlink=True)
#     bpy.data.meshes.remove(pts_mesh, do_unlink=True)
#     bpy.data.materials.remove(mat, do_unlink=True)
#     bpy.data.node_groups.remove(gn, do_unlink=True)


def plotting_in_blender(verts_np, si, mi, color_rgba, frameCount, dot_radius=0.02):
    # --- vertices & faces ---
    verts = [tuple(map(float, v)) for v in np.asarray(verts_np).reshape(-1, 3)]
    faces1 = [tuple(int(v) for v in row) for row in firstfaces]

    # If faces are 1-based, shift to 0-based
    if max(max(f) for f in faces1) == len(verts):
        faces1 = [tuple(v-1 for v in f) for f in faces1]

    if max(max(f) for f in faces1) >= len(verts):
        raise ValueError("Face indices out of range.")

    # --- Material (shared) ---
    rgba = list(color_rgba) if len(color_rgba) == 4 else list(color_rgba) + [1.0]
    mat = bpy.data.materials.new(name=f"CellMat_{mi}_{frameCount}")
    mat.use_nodes = True
    bsdf = mat.node_tree.nodes.get("Principled BSDF") or mat.node_tree.nodes.new("ShaderNodeBsdfPrincipled")
    bsdf.inputs["Base Color"].default_value = rgba
    bsdf.inputs["Alpha"].default_value = 1
    mat.blend_method = 'BLEND'
    # mat.shadow_method = 'HASHED'

    # --- Surface object ---
    surf_mesh = bpy.data.meshes.new(f'cellMesh_{mi}_{frameCount}')
    surf_mesh.from_pydata(verts, [], faces1)
    surf_mesh.update()

    surf_obj = bpy.data.objects.new(f'cellObject_{mi}_{frameCount}', surf_mesh)
    surf_obj.data.materials.append(mat)
    bpy.context.scene.collection.objects.link(surf_obj)

    # --- Use *exact* vertices from surface for the points object ---
    surf_verts = [v.co[:] for v in surf_mesh.vertices]
    pts_mesh = bpy.data.meshes.new(f'pointsMesh_{mi}_{frameCount}')
    pts_mesh.from_pydata(surf_verts, [], [])
    pts_mesh.update()

    pts_obj = bpy.data.objects.new(f'pointsObject_{mi}_{frameCount}', pts_mesh)
    bpy.context.scene.collection.objects.link(pts_obj)

    # --- Sphere to instance ---
    bpy.ops.mesh.primitive_ico_sphere_add(radius=dot_radius, location=(0, 0, 0))
    ico = bpy.context.active_object
    ico.name = f"dot_proto_{mi}_{frameCount}"
    ico.data.materials.append(mat)

    # --- Geometry Nodes: instance sphere on vertices ---
    gn = bpy.data.node_groups.new(f"GN_InstDots_{mi}_{frameCount}", 'GeometryNodeTree')
    mod = pts_obj.modifiers.new(name="GN_InstDots", type='NODES')
    mod.node_group = gn

    nodes = gn.nodes
    links = gn.links
    nodes.clear()

    # Create explicit geometry sockets for Blender 4.5
    gn.interface.new_socket(name="Geometry", in_out='INPUT', socket_type='NodeSocketGeometry')
    gn.interface.new_socket(name="Geometry", in_out='OUTPUT', socket_type='NodeSocketGeometry')

    ng_input = nodes.new("NodeGroupInput")
    ng_output = nodes.new("NodeGroupOutput")

    obj_info = nodes.new('GeometryNodeObjectInfo')
    obj_info.inputs['As Instance'].default_value = True
    obj_info.inputs['Object'].default_value = ico

    inst_on_points = nodes.new('GeometryNodeInstanceOnPoints')

    links.new(ng_input.outputs["Geometry"], inst_on_points.inputs["Points"])
    links.new(obj_info.outputs["Geometry"], inst_on_points.inputs["Instance"])
    links.new(inst_on_points.outputs["Instances"], ng_output.inputs["Geometry"])

    # --- Render ---
    imagename = basename + f"Figure_number_{frameCount}_{si:.2f}.png"
    bpy.context.scene.render.filepath = imagename
    bpy.ops.render.render(write_still=True)

    # --- Cleanup ---
    bpy.data.objects.remove(surf_obj, do_unlink=True)
    bpy.data.objects.remove(pts_obj, do_unlink=True)
    bpy.data.objects.remove(ico, do_unlink=True)
    bpy.data.meshes.remove(surf_mesh, do_unlink=True)
    bpy.data.meshes.remove(pts_mesh, do_unlink=True)
    bpy.data.materials.remove(mat, do_unlink=True)
    bpy.data.node_groups.remove(gn, do_unlink=True)

df_F = pd.read_csv(argv[6][:-4]+".csv")

locate_and_plot_by_shape_index(df_F, argv[6], bins=10, samples=10, nSkip=nSkip)