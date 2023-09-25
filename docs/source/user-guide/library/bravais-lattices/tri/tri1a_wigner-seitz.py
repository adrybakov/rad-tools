import radtools as rad

l = rad.lattice_example("TRI1a")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="wigner-seitz")
# Save an image:
backend.save("tri1a_wigner-seitz.png", elev=9, azim=18, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=9, azim=18)
