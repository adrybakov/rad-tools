import radtools as rad

l = rad.lattice_example("HEX")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="wigner-seitz")
# Save an image:
backend.save("hex_wigner-seitz.png", elev=32, azim=10, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=32, azim=10)
