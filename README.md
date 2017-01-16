# Curved_shape_representation
This README.md file was generated on 2017/01/16/ by Girum G Demisse.

## GENERAL INFORMATION 
This is a MATLAB implementation of curved shape representation
approach presented in:
* Demisse, G.G, Aouada,D, Ottersten, B. *Similarity Metric For Curved
  Shapes in Eculidean Space.*, IEEE CVPR 2016.
* Demisse, G.G, Aouada,D, Ottersten, B. *Deformation Based Curved Shape
  Representation*., Submitted to TPAMI.<br />

## Code of conduct
There is only one requirment for using the software: 
***If you use this software (Curved_shape_representation) in any way or form cite the papers listed above.***
   
## Organization
The software package is organized in three main folders:
* Classes   - includes static and value classes.
* 3rdparty  - includes codes from a third party.
* DP_opt    - includes dynamic programming based optimizations.

## Example dataset
Curved shapes, extracted from [kimia99](http://vision.lems.brown.edu/content/available-software-and-databases), are included with the software package.
* Each coloumn of "kimia99.mat" has 11 elements from the same shape category. In total, there are 9 shape catagories.

## HOW TO:
Three commented demo scripts are included with the package.
* Demo_HOWTO_curve_representation.m - shows how to:
  * Create a curve object-- build representation of a curve.(I recommand you start from here.)
  * Estimate point correspondance, using both unfiorm sampling and optimal sampling, between two Curve objects.
  * Compute geodesic curve and distance between two curved shapes.
  * Plot results.
* Demo_HOWTO_curve_model.m - shows how to:
  * Compute Karcher mean of a set of curved shape representations.
  * Compute K-clusters from a dataset of curved shape representations.
* Demo_HOWTO_deformation_transfer.m - shows how to:
  * Extract deformation that acts from the left, given two curved shape.
  * Transfer a deformation to a given Curved shape.
  
## Author:
*Girum G. Demisse*, girumdemisse@gmail.com

 
