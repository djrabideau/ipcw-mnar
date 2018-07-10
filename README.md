# ipcw-mnar

Description:
--------------
This repository houses code to run simulations for Inverse Probability of Censoring Weighting when data are Missing Not At Random (IPCW-MNAR). [1]  The project was co-authored by Judith J. Lok (Boston University, Boston, MA), Constantin T. Yiannoutsos (Indiana University Fairbanks School of Public Health, Indianapolis, IN), Ronald J. Bosch (Center for Biostatistics in AIDS Research, Harvard T.H. Chan School of Public Health, Boston, MA), and Dustin J. Rabideau (Harvard T.H. Chan School of Public Health, Boston, MA). The code was written by DJR. This public repository was setup in June 2018.

Files:
--------------
batch.txt:	This file was used to run sim-cd4-trajectories.R via command line on a shared cluster that used a Slurm job scheduler (https://slurm.schedmd.com)

simMNAR.R:	This R code contains a function "simMNAR" that simulates MNAR data, generates CD4 trajectories, and compares 4 different analysis methods. This function is called from sim-cd4-trajectories.R.

sim-cd4-trajectories.R:	This R code runs the simulation study in parallel and provides a simple plot of the results. It is primarily setup to be called via the command line, but can be modified to run in the R GUI or RStudio by modifying some stuff up front.

References:
--------------
1. Kiragga AN, Lok JJ, Musick BS, et al. (2014) CD4 trajectory adjusting for dropout among HIV-positive patients receiving combination antiretroviral therapy in an East African HIV care centre. Journal of the International AIDS Society. 17(1).
