# Bloch-McConnell Simulation Code

**MATLAB code for Bloch simulations of FLASH sequences with and without MT-weighting.**

## Requirements:

- Requires the qMRlab toolbox (http://qmrlab.org/) which can be downloaded here: https://github.com/qMRLab/qMRLab
- Can be combined with other sequence simulation code for correcting your maps for B1+ inhomogeneity https://github.com/TardifLab/MTsatB1correction 
- Assumes an MP2RAGE is used for T1 mapping, and thus for the noise calculations. So we need the MP2RAGE code from: https://github.com/JosePMarques/MP2RAGE-related-scripts
- Contains some misc code such as colormaps https://github.com/christopherrowley/NeuroImagingMatlab

***

## Code Overview:

**Not included in the paper, but there is code for RF spoiling, gradient spoiling and spin diffusion here.**

Code is broken up into sections depending on what you are looking to do. 

- **TestCases/** provides some test functions to ensure you are getting expected values
- **Batch_sim/** includes some helper functions for running large simulations and SNR calculation code
- **Functions/** contains the functions for sequence simulations
- **fitTissParams/** has the code and data used to extract the tissue parameters used in the paper.
- **kspaceWeighting/** has the code used for calculating and filling the sampling table, along with scaling based on the MNI atlas.


***


## Instructions/Notes:

- modify `setupIHMTsimPaths.m` to point to the necessary dependencies. 
- call script `TestCases/RunTestCases.m` to make sure the code works, after modifying the save directory for the outputs 
    - Figures are currently already saved for what the results should look like
- Units of time are expressed in seconds, unless stated otherwise.
- If you are going to change the matrix size, you will need to re-reun `kspaceWeighting/Generate_GM_WM_kspaceMatrices.m`

***

## Batch_sim

This is for running simulations over large parameter ranges.  
Use the code in Batch_sim/3T_v1/.  
Change parameters in scripts in 1, and 2, to match what you are interested in. Pay attention to defining the number of nodes used in the cluster on line 20. If unsure, commment out.

1. Run `sim_ihMT_cortex3T_BATCH_conventional.m`
2. Run `sim_ihMT_cortex3T_BATCH_boost_part1.m`  -> simulations
3. Run `sim_ihMT_cortex3T_BATCH_boost_part2.m`  -> signal extractions
4. Run `CR_batch_calculateSNR.m`
5. Run `Batch_sim/functions/FigureGen/plotSimResultsIHMT_figures_SNR_v1.m` 

Depending on how many parameters combinations you set up, this could take a lot of time to generate simulation results.
Expect a few days for perfect spoiling with >200k parameter combos, and a few weeks without perfect spoiling.

***

## Parameter Fitting:
Two options:

1. Run a quick parameter fit using `fitTissParams/SingleFitTesting.m` and you will get a plot showing the fit
2. Run a longer simulation to get better parameters in `fitTissParams/CR_fitTissueParameters.m`

***

Some previously generated files were too large for github. If you would like them, send an email to the address below. The files include boosted simulation results `Boost_dual_Signal.mat`,
`Boost_single_Signal.mat`, and one of the kspaceWeighting matrices `GM_seg_MNI_152_kspace.mat`. The first two will take a while to recreate, the last would be very fast to calculate.


**If used, please reference the following publication:**

Rowley CD, Leppert IR, Campbell JSW, Pike GB, Tardif CL. Optimization of acquisition parameters for cortical inhomogeneous magnetization transfer (ihMT) imaging using a rapid gradient echo readout. Magnetic Resonance in Medicine. 2023 
https://onlinelibrary.wiley.com/doi/10.1002/mrm.29754

**For additional help with this code, contact christopher.rowley@mcgill.ca**

