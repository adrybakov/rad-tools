import radtools as rad

l = rad.lattice_example("TRI1b")
l.plot("primitive")
# Save an image:
l.savefig("tri1b_real.png", elev=12, azim=11, dpi=300, bbox_inches="tight")
# Interactive plot:
l.show(elev=12, azim=11)
