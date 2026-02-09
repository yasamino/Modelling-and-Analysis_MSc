# Hertwig's Rule: Long Axis vs Random Division

This folder contains simulations and analysis comparing the effect of cell division style on the dynamics and progress of the system. We assume two different cell division style: One along the long axis (Hertwig's rule) versus random division orientations.

## Overview

**Long Axis Division (Hertwig's Rule)**

- Division plane is perpendicular to the cell's major axis of elongation.
- Orientation is determined by cell geometry (eigenvalues of gyration tensor).
- Biologically observed in the sarcoma cell line.

**Random Division**

- Division axis is chosen stocastically from all possible orientations in 3D space
- No geometric constraint
- Tests the significance of Hertwig's rule

## Visualizations

### Division Dynamics

<video src="videos/animated_polar_division.mp4" controls width="600"></video>

*Blender render showing division dynamics*

### Gyration Tensor Analysis

#### Major-axis division rule

*Top view (XY)*

<video src="elong_gyration_XY.mp4" controls width="600"></video>

*Side view (YZ)*

<video src="elong_gyration_YZ.mp4" controls width="600"></video>

*show the Gyration tensor analysis of elongation patterns*

#### Random-axis division rule

*Top view (XY)*

<video src="random_gyration_XY.mp4" controls width="600"></video>

### Cell Shape Elongation Analysis

<img src="videos/elong_XY_growth_and_Division_surface_long_axis_4.gif" alt="Vector field visualization of elongation" width="600">

*Shows the vector field visualization of elongation ([Comelles et. al.](https://elifesciences.org/articles/57730)) magnitude and direction*

### Directories

- `data/` - Binary coordinate files
- `scripts/` - Python analysis scripts
- `videos/`
