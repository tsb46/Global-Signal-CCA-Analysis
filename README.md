# Global-Signal-CCA-Analysis
Canonical Correlation Analysis of Global Signal Maps and Behavior

This repository contains the analysis pipeline for a project linking each subject's patterns of Global BOLD signal and behaviora/cognitive measures. In particular, we conduct a canonical correlation analysis linking Global Signal Whole-Brain Beta Maps and subject behavioral measures provided by the Human Connectome Project (N > 1000). Importantly, to follow along with this code, you'll need to download have all the relevant input data: Global Signal Beta Maps and subject behavioral measures (.csv file). Thus, you'll need to run through the appopriate channels with the Human Connectome Project to get that data.

In addition, you'll need the following code dependencies on your path:

Nearest Symmetric Positive Definite Matrix Code:
https://www.mathworks.com/matlabcentral/fileexchange/42885-nearestspd

FSLNETS toolbox:
https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FSLNets

FSL PALM toolbox:
https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/PALM/UserGuide

You will also need some way to load in the cifti files, if your Global Signal Beta Maps are in cifti format (ours were). If you have all of these in .MAT files already than no need for it. We used a matlab cifti toolbox for this:

https://github.com/Washington-University/cifti-matlab
