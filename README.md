# **Description**

C++ programmes to compute thermomechanical properties of solids and liquids from the outputs of LAMMPS molecular dynamics simulations.

**pp_solid.cpp:** Computes how the temperature of a specific section of solid varies in time. Used to ensure equilibration temperature has been attained, or identify if steady state has been achieved.

**pp_liquid.cpp:** Computes how the temperature and pressure of a specific section of the liquid varies in time. Used to ensure equilibration temperature has been attained, or identify if steady state has been achieved.

**pp_bins_solid.cpp:** Divides domain into small one-dimensional subsections (bins) and computes time-averaged temperature in each bin of the solid water. Used to compute temperature profiles, liquid pressures, and near-wall density fluctuations.

**pp_bins_water.cpp:** Divides domain into small one-dimensional subsections (bins) and computes time-averaged temperature, pressure, and density in each bin of water. Used to compute temperature profiles, liquid pressures, and near-wall density fluctuations.

# **Pre-requisites**
LAMMPS dump file, with per-atom information formatted as: 

`x y z vx vy vz c_dstress[1] c_dstress[2] c_dstress[3]`

where dstress is obtained using the LAMMPS command

`compute    stress/atom`

# **Compilation and usage**
In terminal:

`g++ pp_bins_water.cpp`

`./a.out`
