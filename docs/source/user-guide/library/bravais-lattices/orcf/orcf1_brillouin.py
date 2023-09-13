import radtools as rad

l = rad.lattice_example("ORCF1")
l.plot("brillouin-kpath")
# Save an image:
l.savefig("orcf1_brillouin.png", elev=21, azim=49, dpi=300, bbox_inches="tight")
# Interactive plot:
l.show(elev=21, azim=49)