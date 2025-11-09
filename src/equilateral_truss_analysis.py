import numpy as np
import math
import matplotlib.pyplot as plt

# ---------------------------------------------------------------
# USER INPUTS
# ---------------------------------------------------------------
print("Equilateral Triangular Truss Analysis")
try:
    A = float(input("Enter cross-sectional area A (m²): "))
    E = float(input("Enter Young’s Modulus E (Pa): "))
    a = float(input("Enter side length a (m): "))
    P = float(input("Enter load at top node P (N): "))
except ValueError:
    print("Invalid input! Please enter numerical values only.")
    exit()

print("\nSupport Condition for Node 2:")
print("  1 → Fixed support")
print("  2 → Roller support")
node2_type = input("Enter 1 or 2:")

if node2_type not in ["1", "2"]:
    print("Invalid choice! Defaulting to roller support at Node 2.")
    node2_type = "2"

# ---------------------------------------------------------------
# NODE COORDINATES
# ---------------------------------------------------------------
nodes = {
    1: np.array([0.0, 0.0]),
    2: np.array([a, 0.0]),
    3: np.array([a/2.0, math.sqrt(3)/2.0*a])
}

# ELEMENT CONNECTIVITY
elements = {
    1: (1, 2),
    2: (1, 3),
    3: (2, 3)
}

ndof = 2 * len(nodes)
K_global = np.zeros((ndof, ndof))

# ---------------------------------------------------------------
# ELEMENT STIFFNESS MATRIX FUNCTION
# ---------------------------------------------------------------
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
    return k, L, c, s

# ---------------------------------------------------------------
# GLOBAL STIFFNESS MATRIX ASSEMBLY
# ---------------------------------------------------------------
for e, (i, j) in elements.items():
    k, L, c, s = element_stiffness(i, j)
    dof_map = [2*(i-1), 2*(i-1)+1, 2*(j-1), 2*(j-1)+1]
    for a_i in range(4):
        for b_i in range(4):
            K_global[dof_map[a_i], dof_map[b_i]] += k[a_i, b_i]

# ---------------------------------------------------------------
# LOAD VECTOR
# ---------------------------------------------------------------
F = np.zeros(ndof)
F[2*(3-1)+1] = -P  # vertical downward load at node 3

# ---------------------------------------------------------------
# BOUNDARY CONDITIONS
# ---------------------------------------------------------------
fixed_dofs = [0, 1]  # Node 1 fixed

if node2_type == "1":
    fixed_dofs += [2, 3]  # Node 2 fixed
    support_label = "Fixed"
else:
    fixed_dofs += [3]     # Node 2 roller
    support_label = "Roller"

free_dofs = [i for i in range(ndof) if i not in fixed_dofs]

# ---------------------------------------------------------------
# SOLVE FOR DISPLACEMENTS
# ---------------------------------------------------------------
K_ff = K_global[np.ix_(free_dofs, free_dofs)]
K_fr = K_global[np.ix_(free_dofs, fixed_dofs)]
F_f = F[free_dofs]

U_f = np.linalg.solve(K_ff, F_f)
U = np.zeros(ndof)
U[free_dofs] = U_f

# ---------------------------------------------------------------
# REACTIONS
# ---------------------------------------------------------------
R = K_global @ U - F

# ---------------------------------------------------------------
# MEMBER FORCES
# ---------------------------------------------------------------
member_forces = {}
for e, (i, j) in elements.items():
    _, L, c, s = element_stiffness(i, j)
    dof_map = [2*(i-1), 2*(i-1)+1, 2*(j-1), 2*(j-1)+1]
    Ue = U[dof_map]
    T = (E / L) * np.array([-c, -s, c, s])
    F_int = (A * T) @ Ue
    member_forces[e] = F_int

# ---------------------------------------------------------------
# OUTPUT RESULTS
# ---------------------------------------------------------------
print("\n================= RESULTS =================")
print(f"Support Type at Node 2: {support_label}")
print("\nNodal Displacements (in meters):")
for n in range(1, len(nodes)+1):
    print(f"  Node {n}: Ux = {U[2*(n-1)]:.6e}, Uy = {U[2*(n-1)+1]:.6e}")

print("\nReaction Forces (in N):")
for n in [1, 2]:
    print(f"  Node {n}: Rx = {R[2*(n-1)]:.3f}, Ry = {R[2*(n-1)+1]:.3f}")

print("\nMember Axial Forces (in N):")
for e in elements:
    print(f"  Member {e} ({elements[e][0]}-{elements[e][1]}): {member_forces[e]:.3f}")

# ---------------------------------------------------------------
# FBD PLOT
# ---------------------------------------------------------------
plt.figure(figsize=(6,6))
for e, (i, j) in elements.items():
    xi, yi = nodes[i]
    xj, yj = nodes[j]
    plt.plot([xi, xj], [yi, yj], 'k-', lw=2)
    cx, cy = (xi+xj)/2, (yi+yj)/2
    plt.text(cx, cy, f"{e}", color='blue', fontsize=10)

# Supports
plt.plot(nodes[1][0], nodes[1][1], 'rs', label="Fixed Support (Node 1)")
if node2_type == "1":
    plt.plot(nodes[2][0], nodes[2][1], 'rs', label="Fixed Support (Node 2)")
else:
    plt.plot(nodes[2][0], nodes[2][1], 'go', label="Roller Support (Node 2)")

# Load arrow
x3, y3 = nodes[3]
plt.arrow(x3, y3, 0, -0.2*a, head_width=0.05*a, head_length=0.05*a, color='red', label="Applied Load (P)")

plt.text(x3+0.02*a, y3-0.1*a, "P", color='red', fontsize=12)
plt.axis("equal")
plt.grid(True, linestyle="--", alpha=0.5)
plt.legend()
plt.title("Free Body Diagram of Equilateral Triangular Truss")
plt.show()


