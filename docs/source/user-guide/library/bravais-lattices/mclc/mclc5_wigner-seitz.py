import radtools as rad

l = rad.lattice_example("MCLC5")
l.plot("wigner-seitz")
# Save an image:
l.savefig("mclc5_wigner-seitz.png", elev=13, azim=66, dpi=300, bbox_inches="tight")
# Interactive plot:
l.show(elev=13, azim=66)
