import radtools as rad

l = rad.lattice_example("MCLC1")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="brillouin-kpath")
# Save an image:
backend.save("mclc1_brillouin.png", elev=21, azim=53, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=21, azim=53)
