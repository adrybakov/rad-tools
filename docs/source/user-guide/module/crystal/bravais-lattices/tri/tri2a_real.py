import radtools as rad

l = rad.lattice_example("TRI2a")
l.plot("primitive")
# Save an image:
l.savefig("tri2a_real.png", elev=39, azim=44, dpi=300, bbox_inches="tight")
# Interactive plot:
l.show(elev=39, azim=44)
