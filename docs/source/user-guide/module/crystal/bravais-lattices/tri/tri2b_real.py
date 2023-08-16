import radtools as rad

l = rad.lattice_example("TRI2b")
l.plot("primitive")
# Save an image:
l.savefig("tri2b_real.png", elev=17, azim=54, dpi=300, bbox_inches="tight")
# Interactive plot:
l.show(elev=17, azim=54)
