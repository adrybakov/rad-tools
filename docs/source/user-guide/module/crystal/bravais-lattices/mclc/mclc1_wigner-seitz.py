import radtools as rad

l = rad.lattice_example("MCLC1")
l.plot("wigner-seitz")
# Save an image:
l.savefig("mclc1_wigner-seitz.png", elev=45, azim=67, dpi=300, bbox_inches="tight")
# Interactive plot:
l.show(elev=45, azim=67)
