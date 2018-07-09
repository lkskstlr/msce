# Multistage Clonal Expansion (MSCE) MATLAB Simulation Toolbox



## Overview
The repo is structured as follows
+ `src`: All code
 + `msce`: Core functions to simulate MSCE models
 + `tests`: Scripts that are used to generate the figures or test functions
+ `tex`: All tex related files
 + `paper`: Main file explaining the method


## Installation
+ Clone the repository: `git clone https://github.com/lkskstlr/msce.git`
+ Add the src folder to your MATLAB path
+ Install dependencies:
 + **chebfun** is an excellent MATLAB package for Chebyshev polynomial approximation from Oxford's Numerical Analysis Group. It can be obtained from their [website](http://www.chebfun.org/). We recommend installing through git which they call the developer option.
 + **cbrewer** is a MATLAB package that wraps the nice [colorbrewer](https://github.com/altmany/export_fig) colors. It is used for plotting. It can be found on [fileexchange](https://www.mathworks.com/matlabcentral/fileexchange/34087-cbrewer---colorbrewer-schemes-for-matlab).
 + **export_fig** is used for printing figures. The MSCE core should work without it. Find the github repo [here](https://github.com/altmany/export_fig).
