import radtools as rad

l = rad.lattice_example("MCLC3")
l.plot("brillouin-kpath")
# Save an image:
l.savefig("mclc3_brillouin.png", elev=23, azim=34, dpi=300, bbox_inches="tight")
# Interactive plot:
l.show(elev=23, azim=34)
