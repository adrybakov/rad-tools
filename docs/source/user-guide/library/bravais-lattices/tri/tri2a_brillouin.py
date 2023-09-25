import radtools as rad

l = rad.lattice_example("TRI2a")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="brillouin-kpath")
# Save an image:
backend.save("tri2a_brillouin.png", elev=10, azim=32, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=10, azim=32)
