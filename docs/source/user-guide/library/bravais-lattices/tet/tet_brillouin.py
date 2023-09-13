import radtools as rad

l = rad.lattice_example("TET")
l.plot("brillouin-kpath")
# Save an image:
l.savefig("tet_brillouin.png", elev=30, azim=23, dpi=300, bbox_inches="tight")
# Interactive plot:
l.show(elev=30, azim=23)
