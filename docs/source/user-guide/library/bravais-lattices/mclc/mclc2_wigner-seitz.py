import radtools as rad

l = rad.lattice_example("MCLC2")
l.plot("wigner-seitz")
# Save an image:
l.savefig("mclc2_wigner-seitz.png", elev=13, azim=29, dpi=300, bbox_inches="tight")
# Interactive plot:
l.show(elev=13, azim=29)
