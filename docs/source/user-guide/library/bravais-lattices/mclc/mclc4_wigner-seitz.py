import radtools as rad

l = rad.lattice_example("MCLC4")
l.plot("wigner-seitz")
# Save an image:
l.savefig("mclc4_wigner-seitz.png", elev=5, azim=66, dpi=300, bbox_inches="tight")
# Interactive plot:
l.show(elev=5, azim=66)
