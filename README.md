# Modeling and Reverse Time Migration Using Modified Nearly Analytical Discretization
MNAD-MRTM is a MATLAB implementation for acoustic wave modeling and reverse time migration (RTM) using the Modified Nearly Analytical Discretization (MNAD) numerical method.
This repository provides the source code and example scripts for applying MNAD to seismic wave propagation modeling and seismic imaging.

## Features
- Acoustic wave propagation modeling using the MNAD method.
- Reverse time migration (RTM) using the MNAD method.
- Example scripts for seismic modeling and RTM.

## Overview of MNAD Method
A significant challenge in numerical methods arises from inaccuracies in derivatives approximation, making the wave's velocity frequency-dependent and causing numerical dispersion. It is an unphysical artifact that degrades modeling and imaging, particularly at higher frequencies and over time. Avoiding numerical dispersion at high frequencies requires finer spatial grids, which substantially increase computational costs to achieve high-resolution imaging results.
The modified nearly analytical discretization (MNAD) method reduces numerical dispersion by incorporating additional analytical relations through the coupled solution of the wavefield and its gradient fields, using higher-order Taylor expansions, and improving spatial derivative estimation via energy conservation optimization. 
MNAD operates with a compact stencil and maintains stability over long propagation times. Furthermore, it achieves considerable phase-velocity accuracy with only two points per minimum wavelength, requiring a coarser spatial and temporal grid than conventional methods for the same frequency spectrum. These advantages make MNAD particularly beneficial for large-scale studies with long-time recordings and broadband acquisitions.

Table 1 compares the machine-independent computational cost and memory usage of the MNAD with the conventional FD methods. The machine-independent computational cost and memory usage are defined by the total number of operations during the simulation, denoted as N_op, performed at N_x spatial grid points with over N_t time steps, which are controlled by the Courant number and CFL condition:

<img width="662" height="111" alt="table1_forGithub" src="https://github.com/user-attachments/assets/79846e1e-b7e6-47cd-b9bc-b368d014a8c8" />


MNAD is introduced for acoustic RTM in large-scale studies, where compact stencils and coarser grids enable data simulation and high-resolution imaging with significantly reduced computational and memory costs compared to conventional methods. 

## Acknowledgments
This study was funded by the National Science Centre,
Poland, grant no. 2023/51/D/ST10/02603.

## Reproducing the Published Results
This repository contains research code developed for scientific purposes. The implementation is provided as-is and has not been developed as a production seismic-imaging package.
The code supports the numerical methods described in the associated publication.

## Citation
If you use this software in your research, please cite the associated publication:

Zand, T. (2026), Reverse time migration with modified nearly analytical discretization:
Efficient large-scale high-resolution imaging and aliasing suppression, Geophysics, https://doi.org/10.1190/GEO-2025-0665

A CITATION.cff file is included to facilitate citation.

## License
This project is distributed under the license specified in LICENSE.


## Installation
Clone or download the repository and open MATLAB in the
MNAD-MRTM directory. Run:

setup

This adds the required source directories to the MATLAB path.

## Contact
For questions about the implementation, please contact the repository authors.

