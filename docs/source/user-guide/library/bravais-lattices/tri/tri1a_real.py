import radtools as rad

l = rad.lattice_example("TRI1a")
l.plot("primitive")
# Save an image:
l.savefig("tri1a_real.png", elev=31, azim=-20, dpi=300, bbox_inches="tight")
# Interactive plot:
l.show(elev=31, azim=-20)
