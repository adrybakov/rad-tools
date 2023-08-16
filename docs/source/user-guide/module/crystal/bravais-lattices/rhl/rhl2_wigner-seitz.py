import radtools as rad

l = rad.lattice_example("RHL2")
l.plot("wigner-seitz")
# Save an image:
l.savefig("rhl2_wigner-seitz.png", elev=30, azim=-29, dpi=300, bbox_inches="tight")
# Interactive plot:
l.show(elev=30, azim=-29)
