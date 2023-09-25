import radtools as rad

l = rad.lattice_example("RHL1")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="brillouin-kpath")
# Save an image:
backend.save("rhl1_brillouin.png", elev=-41, azim=-13, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=-41, azim=-13)
