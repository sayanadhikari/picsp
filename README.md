[![C/C++ CI](https://github.com/sayanadhikari/picsp/actions/workflows/c-cpp.yml/badge.svg)](https://github.com/sayanadhikari/picsp/actions/workflows/c-cpp.yml)
[![build](https://github.com/sayanadhikari/picsp/actions/workflows/make.yml/badge.svg)](https://github.com/sayanadhikari/picsp/actions/workflows/make.yml)
![PICSP Logo](/images/logo.png)

PICSP - (Particle-in-Cell Simulation of Plasma)
===============================================

*PICSP* is an open source scientific program for simulating plasmas using the *Particle-In-Cell* (PIC) method on a structured mesh. The field quantities are solved using *Direct Solver* and *Gauss-Siedel solver*. The focus is to simulate bounded as well as periodic plasma system.


Contributors
------------

Authors:

- [Sayan Adhikari](https://github.com/sayanadhikari): architecture and data structures, pushers, weighting schemes, overall maintainance.
- [Rakesh Moulick](https://github.com/rakeshmoulick): architecture and data structures, poisson solver, pushers, Collision module.
- [Gunjan Sharma](https://github.com/gunjansharma1019): architecture and data structures, poisson solver, pushers.
- [Rinku Mishra](https://github.com/rinku-mishra): pushers, dust charging implementation.


Installation
------------
#### Prerequisites
1. gcc compiler for C++ (g++)
2. GNU make buildsystem
3. git
4. FFTW3
5. HDF5
5. Python >= 3.0 (Not tested for lower)

#### Procedure
First make a clone of the master branch using the following command
```shell
git clone https://github.com/sayanadhikari207/picsp.git
```
Then enter inside the *PICSP* directory.
```shell
cd picsp
```
Now complile and built the *PICSP* code
```shell
make clean
make all
```
#### Troubleshooting (Installation)
If you have installed **HDF5** library from source, you may get some issues regarding the libraries. Either use *conda* environment or change the path of **HFLAGS** in the makefile.

Use the same method to resolve issues with **FFTW** library.
#### Usage
-----
Upon successful compilation, run the code using following command
```shell
./picsp input.ini
```
Plasma Parameter Setup
----------------------
Edit the input.ini or create your own .ini file using the *input.ini* file format (e.g. *test.ini*) and run the code
```shell
./picsp test.ini
```
