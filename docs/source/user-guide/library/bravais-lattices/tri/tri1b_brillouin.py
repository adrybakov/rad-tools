import radtools as rad

l = rad.lattice_example("TRI1b")
l.plot("brillouin-kpath")
# Save an image:
l.savefig("tri1b_brillouin.png", elev=30, azim=42, dpi=300, bbox_inches="tight")
# Interactive plot:
l.show(elev=30, azim=42)
