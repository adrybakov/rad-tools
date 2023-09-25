import radtools as rad

l = rad.lattice_example("ORCF2")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="brillouin-kpath")
# Save an image:
backend.save("orcf2_brillouin.png", elev=15, azim=36, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=15, azim=36)
