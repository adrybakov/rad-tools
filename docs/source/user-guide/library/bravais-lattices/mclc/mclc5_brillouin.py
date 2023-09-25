import radtools as rad

l = rad.lattice_example("MCLC5")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="brillouin-kpath")
# Save an image:
backend.save("mclc5_brillouin.png", elev=-5, azim=59, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=-5, azim=59)
