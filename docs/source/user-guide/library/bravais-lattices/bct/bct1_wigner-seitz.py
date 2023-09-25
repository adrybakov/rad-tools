import radtools as rad

l = rad.lattice_example("BCT1")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="wigner-seitz")
# Save an image:
backend.save("bct1_wigner-seitz.png", elev=26, azim=59, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=26, azim=59)
