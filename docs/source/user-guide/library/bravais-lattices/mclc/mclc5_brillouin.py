import radtools as rad

l = rad.lattice_example("MCLC5")
l.plot("brillouin-kpath")
# Save an image:
l.savefig("mclc5_brillouin.png", elev=-5, azim=59, dpi=300, bbox_inches="tight")
# Interactive plot:
l.show(elev=-5, azim=59)
