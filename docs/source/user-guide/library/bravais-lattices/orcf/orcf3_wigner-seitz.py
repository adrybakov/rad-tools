import radtools as rad

l = rad.lattice_example("ORCF3")
l.plot("wigner-seitz")
# Save an image:
l.savefig("orcf3_wigner-seitz.png", elev=23, azim=38, dpi=300, bbox_inches="tight")
# Interactive plot:
l.show(elev=23, azim=38)
