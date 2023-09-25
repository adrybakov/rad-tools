import radtools as rad

l = rad.lattice_example("TRI1b")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="wigner-seitz")
# Save an image:
backend.save("tri1b_wigner-seitz.png", elev=21, azim=39, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=21, azim=39)
