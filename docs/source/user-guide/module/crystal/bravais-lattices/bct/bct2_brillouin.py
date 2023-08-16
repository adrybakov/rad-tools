import radtools as rad

l = rad.lattice_example("BCT2")
l.plot("brillouin-kpath")
# Save an image:
l.savefig("bct2_brillouin.png", elev=36, azim=28, dpi=300, bbox_inches="tight")
# Interactive plot:
l.show(elev=36, azim=28)
