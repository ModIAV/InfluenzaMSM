# InfluenzaMSM
This code implements a multiscale model of influenza A virus replication for different MOI conditions developed at the MPI Magdeburg. The current model version is documented in [1]. A detailed description of the intra- and extracellular model is provided in [2] and [3], respectively.

The code provided here allows to simulate the model with different MOI conditions and to conduct the parameter estimation performed in [1]. 

## References
1. Rüdiger D, Kupke SY, Laske T, Zmora P, Reichl U. Multiscale modeling of influenza A virus replication in cell cultures predicts infection dynamics for highly different infection conditions. Unpublished.
2. Heldt FS, Frensing T, Reichl U. Modeling the intracellular dynamics of influenza virus replication to understand the control of viral RNA synthesis. Journal of Virology. 2012;86(15):7806-7817.
3. Heldt FS, Frensing T, Pflugmacher A, Gröpler R, Peschel B, Reichl U. Multiscale modeling of influenza A virus infection supports the development of direct-acting antivirals. PLoS Computational Biology. 2013;9(11):e1003372.

## Requirements
Systems Biology Toolbox for MATLAB by Schmidt and Jirstrand (Bioinformatics, 2006), available at http://www.sbtoolbox2.org/main.php?display=download&menu=download

## Optional programs (for faster simulation and optimization)
- C/C++ compiler: Creates MEX-files for a faster simulation with the SB Toolbox (e.g. MinGW 6.3 C/C++ for Windows or GCC for Linux)

- CVODE solver from SUNDIALS: Simulates MEX-files. Cohen and Hindmarsh (Computers in Physics, 1996), available at https://computation.llnl.gov/projects/sundials/sundials-software

- (f)SSm algorithm: Global optimization of parameters. Egea et al. (Journal of Global Optimization, 2007).

## Running the code and main options
The function `InfluenzaMSM_Main.m` is used for model simulation and parameter estimation. A model simulation can be done by running this script. The following main options are available:
-	Define if a model simulation or parameter estimation (only for MOI 73) should be conducted. Set the variable `p.Task` to either 'simulate' or 'optimize'. 

-	Define the MOI in the variable `p.Ex.Moi`.

-	Choose the optimization algorithm. Set `p.Solver` to 'fminsearch', 'fmincon', or 'fSSm' (fSSm algorithm must be installed).

A simulation of the model with MOI 73, 3 or 10-4 reproduces the main results provided in [1].

## Contributors
The code base was written by Stefan Heldt. Code and model extension was performed by Daniel Rüdiger. 

## Citation
If you use this code, we ask you to cite the appropriate papers in your publication.

