import radtools as rad

l = rad.lattice_example("TRI1b")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="brillouin-kpath")
# Save an image:
backend.save("tri1b_brillouin.png", elev=30, azim=42, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=30, azim=42)
