import radtools as rad

l = rad.lattice_example("TRI1a")
l.plot("brillouin-kpath")
# Save an image:
l.savefig("tri1a_brillouin.png", elev=45, azim=-55, dpi=300, bbox_inches="tight")
# Interactive plot:
l.show(elev=45, azim=-55)
