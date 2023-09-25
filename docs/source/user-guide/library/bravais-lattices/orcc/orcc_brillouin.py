import radtools as rad

l = rad.lattice_example("ORCC")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="brillouin-kpath")
# Save an image:
backend.save("orcc_brillouin.png", elev=22, azim=57, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=22, azim=57)
