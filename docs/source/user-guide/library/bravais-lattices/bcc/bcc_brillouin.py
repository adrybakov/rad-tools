import radtools as rad

l = rad.lattice_example("BCC")
l.plot("brillouin-kpath")
# Save an image:
l.savefig("bcc_brillouin.png", elev=11, azim=25, dpi=300, bbox_inches="tight")
# Interactive plot:
l.show(elev=11, azim=25)
