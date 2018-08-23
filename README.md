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
