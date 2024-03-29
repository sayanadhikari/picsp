[![C/C++ CI](https://github.com/sayanadhikari/picsp/actions/workflows/c-cpp.yml/badge.svg)](https://github.com/sayanadhikari/picsp/actions/workflows/c-cpp.yml)
[![build](https://github.com/sayanadhikari/picsp/actions/workflows/make.yml/badge.svg)](https://github.com/sayanadhikari/picsp/actions/workflows/make.yml)

![PICSP Logo](/images/logo.png)

# PICSP - (Particle-in-Cell Simulation of Plasma)

**PICSP** is an open source scientific program for simulating plasmas using the *Particle-In-Cell* (PIC) method on a 2D structured mesh. The field quantities are solved using *Spectral Solver*, *Direct Solver* and *Gauss-Siedel solver*. The focus is to simulate bounded as well as periodic plasma system.

![configspace_v4_3](https://user-images.githubusercontent.com/39854040/141596731-c973ab56-6bc8-4180-9a06-aa1bb0deb5f5.gif)


## Contributors

Authors:

- [Sayan Adhikari](https://github.com/sayanadhikari), (UiO, Norway): architecture and data structures, pushers, weighting schemes, overall maintainance.
- [Rakesh Moulick](https://github.com/rakeshmoulick), (RC, India): architecture and data structures, poisson solver, pushers, Collision module.
- [Gunjan Sharma](https://github.com/gunjansharma1019), (CPP, India): architecture and data structures, poisson solver.
- [Rupak Mukherjee](https://github.com/RupakMukherjee), (PPPL, USA): parallelization, poisson solver.


## Installation

### Prerequisites
1. [GCC compiler for C++ (g++)](https://gcc.gnu.org/)
2. [GNU make buildsystem](https://www.gnu.org/software/automake/) (Mostly comes with every unix based systems)
3. [git](https://git-scm.com/)
4. [FFTW3](http://www.fftw.org/download.html)
5. [HDF5](https://www.hdfgroup.org/solutions/hdf5/)
5. [Python >= 3.0](https://www.python.org/downloads/) (Not tested for lower)

### Procedure
First make a clone of the ``master`` branch using the following command
```shell
git clone https://github.com/sayanadhikari207/picsp.git
```
Then enter inside the *PICSP* directory.
```shell
cd picsp
```
---
**NOTE**
Ubuntu users can run ``ubuntu.sh`` file to install **FFTW3** and **HDF5** libraries with superuser privilege.
```shell
sudo bash ubuntu.sh
```
---
Now complile and built the *PICSP* code.
```shell
make veryclean
make all
```
### Troubleshooting (Installation)
If you have not installed **HDF5** library from source, you may get some issues regarding the library. Try to use ``apt-get`` for ``Debian`` based machines (e.g. ``Ubuntu``), or any default package manager to install respective libraries. For ``MacOS`` ``Homebrew`` is recommended.

For ``Ubuntu``: 
```shell
sudo apt-get install -y fftw3-dev libhdf5-dev libhdf5-serial-dev
```
For ``MacOs``: 
```shell
brew install fftw hdf5@1.10
```

Use the same method to resolve issues with **FFTW** library.
## Usage

Upon successful compilation, run the code using following command
```shell
./picsp input.ini
```
## Plasma Parameter Setup

Edit the ``input.ini`` or create your own ``.ini`` file using the ``input.ini`` file format (e.g. ``test.ini``) and run the code
```shell
./picsp test.ini
```
## Visualization
To visualize the simulation data, run the python scripts from `scripts` directory.
#### Phase-space Animation
```shell
python3 animatePhaseSpace.py 
```
#### Potential Animation
```shell
python3 animateGridData.py -p phi
```
#### Density Animation
##### Ion
```shell
python3 animateGridData.py -p den.i
```
##### Electron
```shell
python3 animateGridData.py -p den.e
```
#### Time Dependent Data
```shell
python3 timeData.py -p energy
```
