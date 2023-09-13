import radtools as rad

l = rad.lattice_example("ORCC")
l.plot("brillouin-kpath")
# Save an image:
l.savefig("orcc_brillouin.png", elev=22, azim=57, dpi=300, bbox_inches="tight")
# Interactive plot:
l.show(elev=22, azim=57)
