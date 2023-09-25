import radtools as rad

l = rad.lattice_example("BCT2")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="wigner-seitz")
# Save an image:
backend.save("bct2_wigner-seitz.png", elev=41, azim=59, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=41, azim=59)
