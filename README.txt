* INTRODUCTION

The LEE code for simulating Dense Plasma Focus operation under the Lee Model
for GNU ForTran by Jorge A. García Gallardo, CNEA (c) 2019-2020
This code is distributed under the GNU General Public License. For More infomartion read LICENSE.txt

This code is based on Lee's RADPFV6.1b.c over macro version 5.13.6 for Microsoft Excel
with corrections from:

- S. Lee, "Theoretical Basis: Plasma Focus Model (Radiative)-S Lee Model (2008)" http://www.intimal.edu.my/school/fas/UFLF/

- S. Lee and S. H. Saw, "Plasma Focus Numerical Experiments" At the joint ICTP-IAEA Workshop on Dense Magnetized Plasmas and Plasma Diagnostics Trieste 15-26 November 2010

- S. Lee, "Plasma Focus Radiative Model: Review of the Lee Model Code", J. Fusion Energ. (2014) 33:319–335

And validated with the data published in Table 3.10 (pp. 186-187) of

- S. Lee and S. H. Saw, "Plasma Science and Technology for Emerging Economies" Springer Nature (2017)
  "Chapter 3 : The Plasma Focus - Numerical Experiments, Insights and Aplications"

Files included in this distribution:

        README.txt  : This file
        LICENSE.txt : Terms of GNU Public License
        lee.f90     : Fortan 90 source code
        makefile    : Makefile for building binary file
        plot.plt    : GNU plot script for creating plots of the relevant parameters.
        ttf.sh      : bash utility to bin the plots with the phase times
        pf1000/*    : folder with example files for the PF-1000 device
        pf400/*     : folder with example files for the PF-400J device
        nx3/*       : folder with example files for the NX3 device

* NOTES

I. This code has been developed for GNU/LINUX platforms and is expected to compile and run in other UN*X like platforms. Running of this code on Microsoft Windows platforms is not assured.

II. Not all functionalities of the original Lee code are present:

        1. Only H, D or T calculations possible, other elements removed
        2. Tapper capabilities removed.
        3. Variable names were kept when possible
        4. Line labels kept as comments when possible
        5. Original comments kept when possible
        
* SYSTEM REQUIREMENTS

This code requires GNU/fortran and associated libraries to be installed in the system. In addition GNUplot is required for the graphics.

* COMPILATION

In a sh/bash terminal just run

    make

then the binary 'lee' should be created.

* USAGE

First fill a file named <name>.in , for example "data.in", with the parameters of the system to be simulated. These paramenters are the same as in the original Lee code.
Then just run:

    ./lee data.in
  
after a short run, the program will create the following files:

data.param.out      : general information of the simulation
data.out            : main data set (current, voltage, etc) as a function of time
data.temp.out       : temperature data as a function of time
data.energy.out     : energy data as a function of time
data.tics.out       : time of the phases for bining purposes

Then run:

    gnuplot plot.plt

to create the plots of current, temperature, coordinates, etc.



  



