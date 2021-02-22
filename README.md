![PICSP Logo](/images/logo.png)

PICSP - (Particle-in-Cell Simulation of Plasma)
===============================================

*PICSP* is an open source scientific program for simulating plasmas using the *Particle-In-Cell* (PIC) method on a structured mesh. The field quantities are solved using *Direct Solver* and *Gauss-Siedel solver*. The focus is to simulate bounded as well as periodic plasma system.


Contributors
------------

Authors:

- [Sayan Adhikari](mailto:sayan.adhikari@fys.uio.no): architecture and data structures, pushers, weighting schemes, overall maintainance
- Rakesh Moulick (see separate branch): architecture and data structures, poisson solver, Collision module
- [Rinku Mishra](mailto:rinku.mishra@cppipr.res.in): pushers, dust charging implementation
- Gunjan Sharma (see separate branch): Physics Study
- Rupali Paul (see separate branch): Physics Study
- Kishor Deka (see separate branches): Physics Study

Installation
------------
#### Prerequisites
1. gcc compiler for C++ (g++)
2. make buildsystem
3. git

#### Procedure
First make a clone of the master branch using the following command
```shell
git clone https://github.com/sayanadhikari207/picsp.git
```
Then enter inside the *PICSP* directory and verify the contents. 
```shell
cd picsp
ls
```
Then built the  *iniparser*
```shell
cd lib/iniparser
make clean
make
cd ../../
```
Now complile and built the *PICSP* code
```shell
make clean
make
``` 
Usage
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
