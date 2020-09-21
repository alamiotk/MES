# MES

The aim of the task was to simulate heat transfer using the finite element method.
I accomplished the task using the C++ language and the IDE from JetBrains: Clion.

The model assumes a 2D-mesh on the plane, consisting of four-node elements.
There has been used a two-point Gaussian integration scheme. 
The individual points have been described on the basis of Gauss-Legendre quadratures.

To implement the FEM 2D thermal problem there has beed used convection boundary condition.

-----------------------------------------------------------------------------------

ABOUT THE CODE:

The "Node" and "Element" classes are used to create objects that form a mesh.

The "Grid" is my main class where I create a grid and loop over time  and elements. 
I do all the aggregations and calculate the solution to the matrix equation here.

The "UniversalElement" is the class in which I perform all calculations related to integration matrixes H, C, HBC and vector P.

-----------------------------------------------------------------------------------

I have attached a pdf with the entire report.
