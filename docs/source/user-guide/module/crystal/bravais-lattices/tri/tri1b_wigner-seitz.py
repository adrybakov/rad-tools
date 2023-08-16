import radtools as rad

l = rad.lattice_example("TRI1b")
l.plot("wigner-seitz")
# Save an image:
l.savefig("tri1b_wigner-seitz.png", elev=21, azim=39, dpi=300, bbox_inches="tight")
# Interactive plot:
l.show(elev=21, azim=39)
