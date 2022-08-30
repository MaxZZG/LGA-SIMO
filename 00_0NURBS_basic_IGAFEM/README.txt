This is the README file for the Isogeometric Analysis (IGA) Matlab code.
The code supports one, two and three dimensional linear elasticity problems.
Extended IGA for hole, inclusion and crack modelling is also implemented. 
Structural mechanics problems including Euler beam and Kirchhoff plates
(rotation-free formulations) are implemented as well.

The features of the code include:

- Global h-refinement using knot insertion is provided for one and two
dimensional meshes. 

- p- and k-refinement via Geopdes library are also supported.

- Extended IGA which is a Partition of Unity enrichment IGA for 2D traction-free
cracks is also implemented. Level sets are used to detect enriched nodes.
However, enrichment functions are defined in terms of standard geometry i.e.
not in terms of level sets. Holes and inclusions are also implemented.

- An ad hoc implementation for 3D cracks is also provided. In this case,
level sets are used both for enrichment detection and enrichment function evaluation.

- Visualization of displacements, stresses are done in Paraview by exporting the
results to a VTU file. Mesh for visualization purpose is Q4 in 2D and B8 in
3D where B8 denotes tri-linear brick elements.

- Inhomogeneous Dirichlet boundary conditions are treated with the penalty
method, Lagrange multiplier method and the least square method. 

- Numerical integration: 

   (o) elements cut by cracks: sub-triangulation 
   (o) elements cut by circular holes/inclusions: adaptive sub-cells (used in Finite Cell Method)

- Fast assembly of FE stiffness using the triple sparse matrix format
(see files of which names end with FastAssembly in folder iga)

- Tsplines based on Bezier extraction operators. 

- Beams and plates without rotation dofs.

Written by:

 Vinh Phu Nguyen, 
 Delft University of Technology, The Netherlands
 nvinhphu@gmail.com

The code is maintained by the Marie Curie "Initial Training Network - Integrating numerical 
simulation and geometric design technology" (ITN-INSIST)

Contact:
Cosmin Anitescu
Bauhaus University - Weimar, Germany
cosmin.anitescu@uni-weimar.de

Please cite our paper when using this package:

V.P. Nguyen, C. Anitescu, S. Bordas, T. Rabczuk
"Isogeometric Analysis: A review and computer implementation aspects"
Mathematics and Computers in Simulation, submitted

Note: Some of the programs use compiled .mex files, which need to be compiled for your architecture.
The compiled versions for the 64-bit version of OS X, Linux and Windows are included in the C_files directory.
If you would like to compile the files yourself, use the compile.m script on OS X or Linux, or the 
compile_win.m script for Windows.

INSTALLATION

It is recommended to add the igafem substructure to the path by right-clicking on the igafem folder within
MATLAB and selecting "Add to Path" -> "Selected Folders and Subfolders". 

WORKING WITH VTS FILES

To visualize the VTS files produced by some of the scripts, it is recommended to use Paraview. It can be
freely downloaded from: http://www.paraview.org/ To visualize the output in ParaView, use the following
procedure:
    1. Open the .vts or .vtu file in ParaView
    2. Click on the "Apply" button on the left side of the screen
    3. Make sure under "Representation" that "Surface" or "Surface with Edges" is selected
    4. Under "Coloring" select displacement "U" or stresses "sigma" and the approriate measure 
        (i.e. "Magnitude", X-direction displacement, etc.)

To view the distorted shape, one can use the "Wrap by Vector" tool. It is the 9th icon in the 3rd row
of toolbar icons in Paraview 4.3.1 (it looks like a bent rectangle). Then under "Properties" select Vector
"U" and an appropriate "Scale Factor" (usually 10 is a good guess) then click again "Apply". 

There is also a tutorial available, which can be accessed from:
http://www.paraview.org/Wiki/The_ParaView_Tutorial

The package is structured into sub-folders which are described in what follows.

I. Main script files: one file for one problem. The names of the files are already
self explanatory. These files can be found in sub-folder "iga".

*******************************************************************
1. iga1D.m
2. igaLShaped.m
3. igaPlateCircularHole.m
4. igaPlateTension.m
5. igaTBeamLeastSquare.m
6. igaCurvedBeam.m
7. igaCurvedBeamFastAssembly.m
8. igaBracketMultiPatches.m
*******************************************************************

Files that start with "xiga" are for linear elastic fracture mechanics problems. 
These files can be found in folder "xiga".

*******************************************************************
1. xigaEdgeCrack.m
2. xigaEdgeShearCrack.m
3. xigaTwoEdgeCracks.m
4. xigaCenterCrack.m
5. xigaInfiniteCrackModeI.m
6. xigaInfiniteCrackModeII.m
7. xigaInfiniteCrackModeILM.m
8. xigaInfiniteCrackModeILeastSquare.m
9. xigaEdgeCrack3d.m
10. xigaInfiniteCrack3d.m
11. xigaEdgeCrackHalfModel.m

*******************************************************************

Some simple 3D elasticity problems are also included in the "iga" folder:

*******************************************************************
1. iga3dBeam.m
2. igaThickCylinder.m
3. igaPinchedCylinder.m
*******************************************************************

Collocation examples, together with a Galerkin implementation with the same code structure,
are located in the "collocation" folder:

*******************************************************************
1. poisson2d\igaAnnulusGalerkin.m
2. poisson2d\igaAnnulusGrevilleCol.m
3. poisson2d\igaAnnulusSuperCol.m
4. elasticity2d\igaPlateHoleGalerkin.m
5. elasticity2d\igaPlateHoleGalerkinNeu.m
6. elasticity2d\igaPlateHoleGrevilleCol.m
7. elasticity2d\igaPlateHoleGrevilleNeu.m
8. elasticity2d\igaPlateHoleSuperCol.m
9. elasticity2d\igaPlateHoleSuperColNeu.m
10. poisson3d\igaHemisphereGalerkin.m
11. poisson3d\igaHemisphereGrevilleCol.m
12. poisson3d\igaHemisphereSuperCol.m
*******************************************************************

II. Data files: control points, knots, orders. These files can be found in
folder "data".


*******************************************************************
1. plateHoleData.m
2. LShapedData.m
3. LShapedC1Data.m
4. plateC1Data.m
5. plateC2Data.m
*******************************************************************

III. Functions or procedures used in main scripts. They are located in different sub-folders depending on their functionalities.

Most of them have comments and quite self explanatory. I therefore
do not describe them here. For post-processing, the following routines are
used:

*******************************************************************
  1. plotStress1: compute displacements and stresses at nodes of a visualization mesh (Q4 mesh) and export the data to a VTU file.
  2. plotStressXIGA: does just the same thing for XIGA problems.
  3. crackedMesh      : visualize cracked mesh as truly cracked domain (Q4 FEM only)
  4. crackedMeshNURBS : the same but for high-order NURBS
*******************************************************************

IV. C_files, C_files_win: contains C implementation of some NURBS algorithms.

V. FEM problem files: start with "fem". These files can be found in folder
"fem".

These files are FEM codes to solve the same problems solved with IGA.
Used as reference solution to check IGA implementation.
Usually the FEM code reads GMSH mesh files (files with extension *.msh).
The *.msh files are created from the GMSH geometry files *.geo.

VI. Some stand alone scripts used to generate figures in the paper.
These files can be found in folder "examples".

*******************************************************************
1. circle.m:            a full circle by quadratic NURBS curve
2. curves1dWeights.m:   effects of weights on NURBS curves
3. deCasteljauExample.m
4. C_files/NURBSplotter.m: plot B-spline basis  functions 
5. discontinuity1d.m: an example showing the insertion of a knot p+1 times to create a discontinuity 
6. hRefinementExample.m: illustrate h-refinement for 1D B-splines
*******************************************************************

VI. Output files (folder "results")

   Output files include postscript files (*.EPS) and VTU files (*.vtu).
   VTU files can be processed by Paraview to plot displacements, stresses
   contours.

VII. The NURBS toolbox: in folder "nurbs-geopdes". Some routines in this toolbox are used for p- and k-refinement.

VIII. Finite Cell Method: simple implementation of the finite cell method with B-spline basis functions.
 The corresponding files are in folder "finite-cell-method".

IX. Multi-patch IGA code: compatible multi-patch IGA code is resided in folder "@patch2D".
