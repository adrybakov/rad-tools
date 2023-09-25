import radtools as rad

l = rad.lattice_example("RHL1")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="wigner-seitz")
# Save an image:
backend.save("rhl1_wigner-seitz.png", elev=19, azim=-19, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=19, azim=-19)
