"""
Equilateral Triangular Truss Analysis
-------------------------------------
Author : Your Name
Date   : 2025-11-09
"""

import numpy as np
import math
import matplotlib.pyplot as plt

# INPUTS
P = 1.0               # applied load magnitude (units: N)
E = 200e9             # Young's modulus (Pa)
A = 0.01              # cross-sectional area (m^2)
a = 1.0               # side length (m)

# Node coordinates
nodes = {
    1: np.array([0.0, 0.0]),
    2: np.array([a, 0.0]),
    3: np.array([a/2.0, math.sqrt(3)/2.0*a])
}

# Element connectivity
elements = {
    1: (1, 2),
    2: (1, 3),
    3: (2, 3)
}

ndof = 2 * len(nodes)
K_global = np.zeros((ndof, ndof))

def element_stiffness(i, j):
    xi, yi = nodes[i]
    xj, yj = nodes[j]
    L = math.hypot(xj - xi, yj - yi)
    c = (xj - xi) / L
    s = (yj - yi) / L
    k = (E * A / L) * np.array([
        [ c*c,  c*s, -c*c, -c*s],
        [ c*s,  s*s, -c*s, -s*s],
        [-c*c, -c*s,  c*c,  c*s],
        [-c*s, -s*s,  c*s,  s*s]
    ])
    return k

# Assemble global stiffness matrix
for e, (i, j) in elements.items():
    k = element_stiffness(i, j)
    dof_map = [2*(i-1), 2*(i-1)+1, 2*(j-1), 2*(j-1)+1]
    for a_ in range(4):
        for b_ in range(4):
            K_global[dof_map[a_], dof_map[b_]] += k[a_, b_]

# Load vector (vertical load at node 3)
F = np.zeros(ndof)
F[2*(3-1)+1] = -P

# Boundary conditions: node1 fixed (dofs 0,1), node2 v constrained (dof 3)
fixed_dofs = [0, 1, 3]
free_dofs = [i for i in range(ndof) if i not in fixed_dofs]

# Reduce and solve
K_ff = K_global[np.ix_(free_dofs, free_dofs)]
F_f = F[free_dofs]
u_f = np.linalg.solve(K_ff, F_f)

# Full displacement vector
U = np.zeros(ndof)
U[free_dofs] = u_f

# Member forces
member_forces = {}
print("\nMember Forces (positive = tension):")
for e, (i, j) in elements.items():
    xi, yi = nodes[i]
    xj, yj = nodes[j]
    L = math.hypot(xj - xi, yj - yi)
    c = (xj - xi) / L
    s = (yj - yi) / L
    dof_map = [2*(i-1), 2*(i-1)+1, 2*(j-1), 2*(j-1)+1]
    u_e = U[dof_map]
    F_e = (E*A/L) * np.array([-c, -s, c, s]).dot(u_e)
    member_forces[e] = F_e
    nature = "Tension" if F_e > 0 else "Compression"
    print(f"  Member {i}-{j}: {F_e:.6g} N ({nature})")

# Reactions
R = K_global.dot(U) - F
print("\nSupport Reactions:")
print(f"  Node 1: Rx = {R[0]:.6g} N, Ry = {R[1]:.6g} N")
print(f"  Node 2: Rx = {R[2]:.6g} N, Ry = {R[3]:.6g} N")

# Plot free-body diagram
fig, ax = plt.subplots(figsize=(6,6))
ax.set_aspect('equal', 'box')

for e, (i, j) in elements.items():
    xi, yi = nodes[i]
    xj, yj = nodes[j]
    ax.plot([xi, xj], [yi, yj], 'k-', lw=2)
    mx, my = (xi+xj)/2, (yi+yj)/2
    f = member_forces[e]
    text = f"F{i}{j} = {f/abs(P):.3f}·P  ({'T' if f>0 else 'C'})"
    ax.text(mx, my, text, fontsize=9, ha='center', va='center')

for nid, (x,y) in nodes.items():
    ax.plot(x, y, 'ro')
    ax.text(x+0.02, y-0.05, f"Node {nid}", fontsize=10)

ax.text(nodes[1][0]-0.1, nodes[1][1]-0.1, "Fixed", color='blue')
ax.text(nodes[2][0]+0.05, nodes[2][1]-0.1, "Roller", color='blue')

x3, y3 = nodes[3]
ax.arrow(x3, y3, 0, -0.2, head_width=0.03, color='red')
ax.text(x3-0.05, y3-0.25, "P ↓", color='red')

ax.axis('off')
ax.set_title("Free-Body Diagram of Equilateral Truss")
plt.show()
