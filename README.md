# AirGun1D

This repository contains MATLAB code for simulating a 1D air gun coupled to an oscillating bubble. The model calculates fields inside the air gun chamber and the acoustic pressure in the water.

This code was used in the paper of Watson, Werpers and Dunham (2019) What controls the initial peak of an air gun source signature? *Geophysics*, 84(2), P27-P45. doi: [10.1190/geo2018-0298.1](https://doi.org/10.1190/geo2018-0298.1) and this repository contains the code necessary to reproduce the figures in this paper.

### How do I get set up? ###
* Clone this repository to your local directory.
* The 1D air gun code is a finite-difference code that uses summation-by-parts (SBP) operators. In order to run the code you will need to link to also clone the [SBPLib repository](https://bitbucket.org/sbpteam/sbplib) and add that to your MATLAB path.
* Once SBPLib is added to your MATLAB path you should be able to run the code and start generating figures.
* There is an example script file in the main directory that will simulate air gun dynamics and plot the bubble radius, acoustic pressure, and pressure inside the source (Figure 3 in the User Guide).
* Script files for generating the figures of Watson, Werpers and Dunham (2018) are contained in the MakeFigs directory. These script files should provide guidance for how to run the 1D air gun code and display outputs.

### Organization of the repository ###
The repository contains several folders:
* MakeFigs - contains the script files for generating each figure.
* SBPSAT - 1D air gun code
* SeismicAirgunCode - lumped parameter model
* sbplib - summation-by-parts operators that are required to run code. 
* Data - data from Lake Seneca field tests. For more details see [Chelminski et al. (2016)](https://www.epmag.com/low-pressure-source-840586#p=full), [Ronen and Chelminski (2018)](http://earthdoc.eage.org/publication/publicationdetails/?publication=92131) or [Ronen et al. (2019)](https://onlinelibrary.wiley.com/doi/abs/10.1111/1365-2478.12774).
* BubbleCode - bubble models that are compared between in Appendix 1.

### 1D Air Gun Code ###
The 1D air gun code is contained in the directory SBPSAT. For details on the code please refer to the paper and in particular Appendix B, which describes the numerical scheme. The main function files are:
* runEulerCode.m - loads the discretization of the air gun and bubble properties and then solves the resulting set of ordinary differential equations using ode45
* DiscrAirgun.m - defines the discretization of the air gun. This function applies all of the boundary treatmens and discretization operators. 
* configAirgun.m - air gun, bubble and ambient properties are defined in this file. For example, the ambient density, gas constant and firing duration are defined in this file. Other parameters, like air gun depth and operating pressure, are defined in the script files.
* bubbleRHS.m - describes the bubble governing equations. Parameters unique to the bubble, like the surface magnification factor, are defined here.

### SBPLIB (Summation-by-Parts Library) ###
* sbplib was developed by Jonatan Werpers, Martin Almquist, Ylva Rydin and Vider Stiernström. The main sbplib repository is [SBPLib repository](https://bitbucket.org/sbpteam/sbplib). Note that there are compatibility issues between the latest version of sbplib and the sbp operators used by AirGun1D. Therefore, it is recommended to use the version of sbplib that is included in the AirGun1D repository.

### Who do I talk to? ###
* Leighton Watson: leightonwatson@stanford.edu, lwatson2@uoregon.edu
* Jonatan Werpers: jonatan.werpers@it.uu.se
