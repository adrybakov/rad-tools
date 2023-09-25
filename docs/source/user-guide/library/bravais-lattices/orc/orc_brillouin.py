import radtools as rad

l = rad.lattice_example("ORC")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="brillouin-kpath")
# Save an image:
backend.save("orc_brillouin.png", elev=35, azim=34, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=35, azim=34)
