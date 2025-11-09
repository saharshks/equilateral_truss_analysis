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

# How to run the python code:-
  Download the Project

Click the green “Code” button.

Choose “Download ZIP”.

After it downloads, extract (unzip) the file on your computer — for example, on your Desktop.

 Open the Folder

Open the folder you just extracted.

Inside, you’ll see your Python file (for example, equilateral_truss_analysis.py).

 Install Python

Make sure you have Python 3 installed.
You can check by opening Command Prompt and typing:

python --version


If you don’t have it, download it.
.

 Install Required Packages

In Command Prompt, go to your project folder using cd. Example:

cd Desktop/equilateral_truss_analysis-main


Then type:

pip install numpy matplotlib

 Run the Code

Now run the program using:

python equilateral_truss_analysis.py

 Enter the Inputs

When the program runs, it will ask for:

Area (A)

Young’s Modulus (E)

Side length (a)

Load at the top node (P)

Type of support at Node 2 (1 = Fixed, 2 = Roller)

After entering these, the program will:

Calculate displacements, member forces, and reactions

Show a Free Body Diagram of the truss
