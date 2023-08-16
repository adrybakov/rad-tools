import radtools as rad

l = rad.lattice_example("ORCI")
l.plot("wigner-seitz")
# Save an image:
l.savefig("orci_wigner-seitz.png", elev=30, azim=12, dpi=300, bbox_inches="tight")
# Interactive plot:
l.show(elev=30, azim=12)
