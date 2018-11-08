# AirGun1D

This repository contains MATLAB code for simulating a 1D air gun coupled to an oscillating bubble. The model calculates fields inside the air gun chamber and the acoustic pressure in the water.

This code was used in the paper of Watson, Werpers and Dunham (2018) What controls the initial peak of an air gun source signature? *Geophysics* and this repository contains the code necessary to reproduce the figures in this paper.


### How do I get set up? ###
* Clone this repository to your local directory.
* The 1D air gun code is a finite-difference code that uses summation-by-parts (SBP) operators. Therefore, in order to run the code you will need to link to also clone the [SBPLib repository](https://bitbucket.org/sbpteam/sbplib) and add that to your MATLAB path.
* Once SBPLib is added to your MATLAB path you should be able to run the code and start generating figures.
* Script files for generating the figures of Watson, Werpers and Dunham (2018) are contained in the MakeFigs directory. This script files should providev guidance for how to run the 1D air gun code.

### Organization of the repository ###
The repository contains several folders:
* MakeFigs - contains the script files for generating each figure.
* SBPSAT - 1D air gun code
* SeismicAirgunCode - lumped parameter model
* Data - data from Lake Seneca field tests. For more details see [Chelminski et al. (2016)](https://www.epmag.com/low-pressure-source-840586#p=full) or [Ronen and Chelminski (2018)](http://earthdoc.eage.org/publication/publicationdetails/?publication=92131).
* BubbleCode - bubble models that are compared between in Appendix 1.

### Who do I talk to? ###

* Leighton Watson: leightonwatson@stanford.edu
* Jonatan Werpers: jonatan.werpers@it.uu.se
