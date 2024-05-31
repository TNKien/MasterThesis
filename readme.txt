READ ME

This folder contains all the code used for the master thesis : 
"Brain Cnnectivity Investigation During Seizures in the ICU"

1) Folder "Data crises" : contains the EEG time series recorded on each of the 8 patients

2) Folder "f_PlotEEG_BrainNetwork" : contains all the necessary functions to plot the brain network and the local clustering coefficient.

3) Folder "main" : contains the main.m file and the necessary functions to build connectivity matrices and calculate connectivity metrics from EEG data

4) Folder "measures" : contains the connectivity measures for the computations of the network edgess

5) Folder "BCT" : contains a wide list of connectivity metrics

6) Folder "Results" : contains the connectivity matrices for each patient, frequency bands, connectivity measure used, time point within the one-hour window AND the corresponding connectivity metrics

7) EEG_plot.m : plot the EEG time series recorded by each channels (21)

8) Lineplot_perPatient.m : plot the temporal evolution of the four connectivity metrics within the one-hour window

9) NP.m : plot the brain network during an interictal and an ictal period

10) stat_test.m : function performing the permutation test on 2 data sets