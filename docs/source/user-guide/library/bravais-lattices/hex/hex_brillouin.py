import radtools as rad

l = rad.lattice_example("HEX")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="brillouin-kpath")
# Save an image:
backend.save("hex_brillouin.png", elev=19, azim=20, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=19, azim=20)
