# Folding pointwise

## The method

The function fs_voxel_features.m contains the MATLAB script of the method described in 
[Local Morphological Measures Confirm that Folding within Small Partitions of the Human Cortex Follows Universal Scaling Law](https://arxiv.org/abs/2103.14061).
For a given subject, it computes the variables average cortical thickness, total surface area, 
exposed surface area, and integrated Gaussian curvature in small regions around each point on the cortex. 
It corrects the surface area measures using the curvature and then computes the new morphological variables K, I, and S.
The method enables the near-pointwise analysis of cortical folding.

The current version excludes points located on the midline or in the insula.

## Input and dependencies

The input of the function is the path to a subject's FreeSurfer folder. The method requires the FreeSurfer files ?h.pial, ?h.thickness, and ?h.pial-outer-smoothed 
(see the [FreeSurfer lGI website](https://surfer.nmr.mgh.harvard.edu/fswiki/LGI)).
The output is saved in the 'surf' folder.

It also needs the MATLAB libraries [CorticalFoldingAnalysisTools](https://github.com/cnnp-lab/CorticalFoldingAnalysisTools) and [iso2mesh](http://iso2mesh.sourceforge.net).