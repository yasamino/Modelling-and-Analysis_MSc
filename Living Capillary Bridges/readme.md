# Living Capillary Bridges  

**Project repository for code, analysis, and simulations**  
*Related paper:* https://doi.org/10.48550/arXiv.2510.14518

**Authors (as listed on the paper):**  
Tytti Kärki, Senna Luntama, Yasamin Modabber, Saila Pönkä, Gonca Erdemci-Tandogan, Mikko Karttunen, Grégory Beaune, Jaakko V. I. Timonen

---

## Overview

This repository contains the code, data, and analysis scripts used in the master's project supporting the paper *Living Capillary Bridges*.  

The goal of this work is to **model, simulate, and analyze the behavior of living capillary bridges** using computational methods. We use [CellSim3D](https://github.com/SoftSimu/CellSim3D), a coarse-grained software and model to simulate cellular systems in three dimensions. The scripts implement key parts of the analysis and visualizations.

---

## Folder Structure
```
├── data/ # data files
├── videos/ # Simulation videos
├── Hertwig's rule/ # Scripts, videos, and discussion regarding the division rules
├── Population analysis/ # Analysis scripts
├── Stochasticity
└── README.md
```
---

## Simulation Videos

The `videos/` directory contains visualisations illustrating the formation and stability of living capillary bridges under different initial conditions.

### Large, stable capillary bridge (R = 100)

- **Side view:** [R100_Sideview.mp4](videos/R100_Sideview.mp4)  
- **Top view:** [R100_Topview.mp4](videos/R100_Topview.mp4)  

These simulations show a large, stable capillary bridge maintaining structural integrity over time.

---

### Small, unstable capillary bridge (R = 70)

- **Side view:** [R70_Sideview.mp4](videos/R70_Sideview.mp4)  
- **Top view:** [R70_Topview.mp4](videos/R70_Topview.mp4)  

These simulations demonstrate an unstable capillary bridge undergoing deformation and eventual loss of stability.
