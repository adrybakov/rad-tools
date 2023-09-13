import radtools as rad

l = rad.lattice_example("RHL2")
l.plot("brillouin-kpath")
# Save an image:
l.savefig("rhl2_brillouin.png", elev=14, azim=-85, dpi=300, bbox_inches="tight")
# Interactive plot:
l.show(elev=14, azim=-85)
