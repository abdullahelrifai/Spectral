# **Description**

C++ programmes to compute vibrational properties of solids and liquids from the outputs of LAMMPS molecular dynamics simulations.

**VDOS.cpp:** Computes the vibrational density of states (VDOS) based on the velocity auto-correlation function of a group of atoms. Is capable of identifying atoms that remain within a specific domain for the duration of the correlation, thus able to compute the VDOS of liquids.

**SDHF.cpp:** Computes the spectral decomposition of heat flux (SDHF) based on the force-velocity cross-correlation function of a group of atoms.

# **Pre-requisites**
Both programmes compute Fourier transforms of time-series data using the [FFTW3 library](https://www.fftw.org/). This must be installed to be included in the compilation of the code.

The formatting of the LAMMPS dump files for each programme is detailed upon below:
### VDOS.cpp
LAMMPS dump file of the group of atoms of interest, with per-atom information formatted as: 

`x y z vx vy vz`

### SDHF.cpp
LAMMPS dump file obtained via the LAMMPS `rerun` command, with per-atom information formatted as: 

`vx vy vz fx fy fz`

# **Compilation and usage**
In terminal:

`g++ <path to FFTW3 library> SDHF.cpp`

`./a.out`

# Validations
### VDOS.cpp
![image](https://github.com/abdullahelrifai/Spectral/assets/160526058/c3551e61-3270-4a09-8de2-7f088cd95c1e)

### SDHF.cpp
![image](https://github.com/abdullahelrifai/Spectral/assets/160526058/949b0707-c2c1-4fb2-9601-0dea80c5f5ed)

