# Curved_shape_representation
This README.md file was generated on 2017/01/16/ by Girum G Demisse.

## GENERAL INFORMATION 
The software is the implementation of curved shape representation
and modelling discussed in:
* Demisse, G.G, Aouada,D, Ottersten, B. *Similarity Metric For Curved
  Shapes in Eculidean Space.*, IEEE CVPR 2016.
* Demisse, G.G, Aouada,D, Ottersten, B. *Deformation Based Curved Shape
  Representation*. Submitted to TPAMI.<br />
*If you use the software in any way or form cite the papers listed above.*
   
## Organization
There are three main folders:
* Classes  - includes static and value classes.
* 3rdparty - includes third party codes.
* DP_opt   - includes dynamic programming based optimizations.

## Example dataset
Curved shapes, extracted from kimia99, are included.
* Each coloumn of "kimia99.mat" has 11 elements from the same category.

## HOW TO:
Three commented demo scripts are included with the package.
* Demo_HOWTO_curve_representation.m - shows how to:
  * Create a curve object-- build representation of a curve.(I recommand you to start from here.)
  * Estimate point correspondance, using both unfiorm sampling and optimal sampling, between two Curve objects.
  * Compute geodesic curve and distance between two curved shapes.
  * Plot results.
* Demo_HOWTO_curve_model.m - shows how to:
  * Compute Karcher mean of a set of curved shape representations.
  * Compute K-clusters from a dataset of curved shape representations.
* Demo_HOWTO_deformation_transfer.m - shows how to:
  * Extract deformation that acts from the left, given two curved shape.
  * Transfer a deformation to a given Curved shape.
## Contributor/s:
Currently, there is only one contributor and main contact: *Girum G. Demisse*, girumdemisse@gmail.com

## License 
 
