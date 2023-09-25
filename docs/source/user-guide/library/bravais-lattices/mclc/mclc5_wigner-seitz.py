import radtools as rad

l = rad.lattice_example("MCLC5")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="wigner-seitz")
# Save an image:
backend.save("mclc5_wigner-seitz.png", elev=13, azim=66, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=13, azim=66)
