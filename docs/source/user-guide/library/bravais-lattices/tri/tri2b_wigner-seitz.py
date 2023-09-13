import radtools as rad

l = rad.lattice_example("TRI2b")
l.plot("wigner-seitz")
# Save an image:
l.savefig("tri2b_wigner-seitz.png", elev=19, azim=13, dpi=300, bbox_inches="tight")
# Interactive plot:
l.show(elev=19, azim=13)
