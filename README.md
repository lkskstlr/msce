# Colorectal Cancer (CRC) MATLAB Simulation Toolbox
---
This toolbox implements the paper
```
Jihyoun Jeon et al. “Evaluation of screening strategies for pre-malignant lesions using a biomathematical approach”. In: Mathematical biosciences 213.1 (2008), pp. 56–70.
```
which can be found on [google scholar](http://scholar.google.com/scholar?btnG=Search%2BScholar&as_q=%22Evaluation%2Bof%2Bscreening%2Bstrategies%2Bfor%2Bpre-malignant%2Blesions%2Busing%2Ba%2Bbiomathematical%2Bapproach.%22&as_sauthors=Jeon&as_occt=any&as_epq=&as_oq=&as_eq=&as_publication=&as_ylo=&as_yhi=&as_sdtAAP=1&as_sdtp=1). The file `tex/paper/paper.pdf` explains an approximation used which decreases the computational time by two to three orders of magnitude while retaining comparable accuracy.

## Overview
The repo is structured as follows
+ `src`: All code
+ `tex`: All tex related files


## Installation
+ Clone the repository: `git clone https://github.com/lkskstlr/msce.git`
+ Add the src folder to your MATLAB path
+ Install dependencies:
  + **chebfun** is an excellent MATLAB package for Chebyshev polynomial approximation from Oxford's Numerical Analysis Group. It can be obtained from their [website](http://www.chebfun.org/). We recommend installing through git which they call the developer option.
  + **cbrewer** is a MATLAB package that wraps the nice [colorbrewer](https://github.com/altmany/export_fig) colors. It is used for plotting. It can be found on [fileexchange](https://www.mathworks.com/matlabcentral/fileexchange/34087-cbrewer---colorbrewer-schemes-for-matlab).
  + **export_fig** is used for printing figures. The MSCE core should work without it. Find the github repo [here](https://github.com/altmany/export_fig).

## Code Structure
+ `src/tests` contains the MATLAB unit tests which can *and should* be run by running `src/tests/crc_run_all_tests.m`.
+ The prefix `crc` indicates functions specific for the Colorectal Cancer paper.
+ Function Nomenclature:
  + `crc_survival*.m` are the survival functions for the `*`-stage model. They return the value.
  + `crc_mksurvival*.m` are the survival functions for the `*`-stage model. They return a chebfun approximation (accurate to nearly machine precision: 1e-15) of the survival function. `crc_mklogsurvival*.m` returns a chebfun for the log of the survival function.`crc_mksurvival*m1.m` returns a chebfun for the survival function minus 1.
  + The above holds also for the hazard functions.
