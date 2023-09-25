import radtools as rad

l = rad.lattice_example("ORC")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="wigner-seitz")
# Save an image:
backend.save("orc_wigner-seitz.png", elev=20, azim=30, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=20, azim=30)
