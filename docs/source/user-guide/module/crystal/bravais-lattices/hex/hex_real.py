import radtools as rad

l = rad.lattice_example("HEX")
l.plot("primitive")
# Save an image:
l.savefig("hex_real.png", elev=35, azim=23, dpi=300, bbox_inches="tight")
# Interactive plot:
l.show(elev=35, azim=23)
