import radtools as rad

l = rad.lattice_example("TRI2b")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="brillouin-kpath")
# Save an image:
backend.save("tri2b_brillouin.png", elev=40, azim=1, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=40, azim=1)
