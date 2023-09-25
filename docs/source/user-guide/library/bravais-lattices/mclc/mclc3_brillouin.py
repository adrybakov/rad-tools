import radtools as rad

l = rad.lattice_example("MCLC3")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="brillouin-kpath")
# Save an image:
backend.save("mclc3_brillouin.png", elev=23, azim=34, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=23, azim=34)
