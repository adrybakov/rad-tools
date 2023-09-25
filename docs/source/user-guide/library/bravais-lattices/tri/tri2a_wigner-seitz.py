import radtools as rad

l = rad.lattice_example("TRI2a")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="wigner-seitz")
# Save an image:
backend.save("tri2a_wigner-seitz.png", elev=30, azim=62, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=30, azim=62)
