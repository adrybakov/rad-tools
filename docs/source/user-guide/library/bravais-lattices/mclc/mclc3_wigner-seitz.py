import radtools as rad

l = rad.lattice_example("MCLC3")
l.plot("wigner-seitz")
# Save an image:
l.savefig("mclc3_wigner-seitz.png", elev=41, azim=61, dpi=300, bbox_inches="tight")
# Interactive plot:
l.show(elev=41, azim=61)
