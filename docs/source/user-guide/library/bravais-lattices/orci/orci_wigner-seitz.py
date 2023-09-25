import radtools as rad

l = rad.lattice_example("ORCI")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="wigner-seitz")
# Save an image:
backend.save("orci_wigner-seitz.png", elev=30, azim=12, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=30, azim=12)
