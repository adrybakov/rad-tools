import radtools as rad

l = rad.lattice_example("TET")
l.plot("wigner-seitz")
# Save an image:
l.savefig("tet_wigner-seitz.png", elev=30, azim=30, dpi=300, bbox_inches="tight")
# Interactive plot:
l.show(elev=30, azim=30)
