import radtools as rad

l = rad.lattice_example("MCL")
l.plot("wigner-seitz")
# Save an image:
l.savefig("mcl_wigner-seitz.png", elev=11, azim=37, dpi=300, bbox_inches="tight")
# Interactive plot:
l.show(elev=11, azim=37)
