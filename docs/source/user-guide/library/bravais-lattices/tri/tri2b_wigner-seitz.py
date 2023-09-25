import radtools as rad

l = rad.lattice_example("TRI2b")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="wigner-seitz")
# Save an image:
backend.save("tri2b_wigner-seitz.png", elev=19, azim=13, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=19, azim=13)
