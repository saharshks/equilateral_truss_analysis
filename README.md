# This Python program performs structural analysis of an equilateral triangular truss
It calculates nodal displacements, member axial forces, and support reactions, and visualizes the structure with a Free Body Diagram

An equilateral truss consists of three members forming a triangle.

The base is supported:

Left node (Node 1) → Always fixed

Right node (Node 2) → Can be roller or fixed, based on user choice

The top node (Node 3) is loaded vertically downward with a given load P.

The code builds the global stiffness matrix, applies boundary conditions, and solves the system:

[K]{U}={F}

to determine all displacements, reactions, and internal member forces.
