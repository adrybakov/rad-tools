import radtools as rad

l = rad.lattice_example("ORCF1")
l.plot("wigner-seitz")
# Save an image:
l.savefig("orcf1_wigner-seitz.png", elev=44, azim=28, dpi=300, bbox_inches="tight")
# Interactive plot:
l.show(elev=44, azim=28)
