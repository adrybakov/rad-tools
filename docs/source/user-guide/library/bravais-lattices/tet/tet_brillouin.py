import radtools as rad

l = rad.lattice_example("TET")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="brillouin-kpath")
# Save an image:
backend.save("tet_brillouin.png", elev=30, azim=23, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=30, azim=23)
