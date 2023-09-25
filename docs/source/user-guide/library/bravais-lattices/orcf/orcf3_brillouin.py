import radtools as rad

l = rad.lattice_example("ORCF3")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="brillouin-kpath")
# Save an image:
backend.save("orcf3_brillouin.png", elev=25, azim=62, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=25, azim=62)
