import radtools as rad

l = rad.lattice_example("MCLC2")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="brillouin-kpath")
# Save an image:
backend.save("mclc2_brillouin.png", elev=11, azim=64, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=11, azim=64)
