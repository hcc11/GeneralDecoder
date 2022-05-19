This folder contains Matlab (R2020a) and C codes for the manuscript: 
A.M. Ni, C. Huang, B. Doiron and M.R. Cohen (2022) A general decoding strategy explains therelationship between behavior and correlated variability. eLife. (in press)

++++++++++++++++++++++++++++++++

To use, first compile all the C codes with the mex compiler in Matlab. Specifically, in the Matlab command line, run the following commands:
mex EIF1DRFfastslowSyn.c
mex spktime2count.c

Sim_Ori_gabor_L2.m simulates a three-layer network with different attentional modulation. Simulation code for Fig. 2.  Saves spike counts.   

CollectSpk.m collect spike counts from every 500 simulations into one file, which will be used to compute model statistics, Fisher information and decoder performance.  

ModelStat.m  Compute model statistics for different attentional modulation (Fig 2D-F)  

FIdecoder_cluster_L1.m computes Fisher information for the specific decoder (Fig. 3B). Calls function FIdecoder.m 

GD_th_cluster.m computes Fisher information for the general decoder (Fig. 3B). 

CollectData.m   Collects Fisher information vs. number of neurons data from different parameter sets.

Local_global_svm.m   Compare performance of specific and general decoders with small number of neurons (Fig. 3D)

RF2D3layer.m is the main simulation function. It contains default parameter values and uses the mex file EIF1DRFfastslowSyn.c for integration. 

genXspk.m generate spike trains of input from Layer 1. 

gen_weights.m generates weight matrices without tuning-dependent connections 

ori_map.m generates a columnar orientation map. 

spktime2count.c converts spike time data to spike counts. 

raster2D_ani.m generates movie of spike rasters from 2D spatial networks. 


+++++++++++++++++++++++++++++++++++

One simulation of a three-layer network for 20 sec takes about 2 hours CPU time and under 3 gb memory. 
(Sim_Ori_gabor_L2.m)



