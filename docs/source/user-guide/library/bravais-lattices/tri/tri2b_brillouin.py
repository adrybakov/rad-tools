import radtools as rad

l = rad.lattice_example("TRI2b")
l.plot("brillouin-kpath")
# Save an image:
l.savefig("tri2b_brillouin.png", elev=40, azim=1, dpi=300, bbox_inches="tight")
# Interactive plot:
l.show(elev=40, azim=1)
