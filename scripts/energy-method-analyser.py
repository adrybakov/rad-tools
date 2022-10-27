atom_template = {
    "J1": [(1, 0), (1, 1), (0, 1), (-1, 0), (-1, -1), (0, -1)],
    "J2": [(2, 1), (1, 2), (-1, 1), (-2, -1), (-1, -2), (1, -1)],
    "J3": [(2, 0), (2, 2), (0, 2), (-2, 0), (-2, -2), (0, -2)]
}

coefficients = {}

for key in atom_template:
    coefficients[key] = 0

supercell = [
    (0, 0),
    (0, 1),
    (0, 2),
    (1, 0),
    (1, 1),
    (1, 2),
    (2, 0),
    (2, 1),
    (2, 2),
]

magnetic_structure = list(map(int, input().split()))


for i_a, atom1 in enumerate(supercell):
    atom1_spin = magnetic_structure[i_a]
    for parameter in atom_template:
        for j_a, atom2 in enumerate(atom_template[parameter]):
            R_i = (atom1[0] + atom2[0]) % 3
            R_j = (atom1[1] + atom2[1]) % 3
            atom2_spin = magnetic_structure[supercell.index((R_i, R_j))]
            coefficients[parameter] += atom1_spin * atom2_spin

print("E = E0 ", end="")
for key, value in coefficients:
    if value >= 0:
        print(f"+ {abs(value)} * {key} ", end="")
    else:
        print(f"- {abs(value)} * {key} ", end="")
