import radtools as rad

l = rad.lattice_example("ORCI")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="brillouin-kpath")
# Save an image:
backend.save("orci_brillouin.png", elev=35, azim=23, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=35, azim=23)
