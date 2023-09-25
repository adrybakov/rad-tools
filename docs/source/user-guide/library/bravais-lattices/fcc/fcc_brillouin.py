import radtools as rad

l = rad.lattice_example("FCC")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="brillouin-kpath")
# Save an image:
backend.save("fcc_brillouin.png", elev=23, azim=28, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=23, azim=28)
