import radtools as rad

l = rad.lattice_example("ORCI")
l.plot("brillouin-kpath")
# Save an image:
l.savefig("orci_brillouin.png", elev=35, azim=23, dpi=300, bbox_inches="tight")
# Interactive plot:
l.show(elev=35, azim=23)
