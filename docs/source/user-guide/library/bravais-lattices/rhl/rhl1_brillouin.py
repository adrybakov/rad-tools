import radtools as rad

l = rad.lattice_example("RHL1")
l.plot("brillouin-kpath")
# Save an image:
l.savefig("rhl1_brillouin.png", elev=-41, azim=-13, dpi=300, bbox_inches="tight")
# Interactive plot:
l.show(elev=-41, azim=-13)
