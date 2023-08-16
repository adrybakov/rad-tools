import radtools as rad

l = rad.lattice_example("ORCF2")
l.plot("wigner-seitz")
# Save an image:
l.savefig("orcf2_wigner-seitz.png", elev=38, azim=14, dpi=300, bbox_inches="tight")
# Interactive plot:
l.show(elev=38, azim=14)
