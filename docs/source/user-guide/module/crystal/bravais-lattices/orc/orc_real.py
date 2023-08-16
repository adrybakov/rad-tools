import radtools as rad

l = rad.lattice_example("ORC")
l.plot("primitive")
# Save an image:
l.savefig("orc_real.png", elev=36, azim=35, dpi=300, bbox_inches="tight")
# Interactive plot:
l.show(elev=36, azim=35)
