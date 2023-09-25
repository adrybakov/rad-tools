import radtools as rad

l = rad.lattice_example("MCLC2")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="wigner-seitz")
# Save an image:
backend.save("mclc2_wigner-seitz.png", elev=13, azim=29, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=13, azim=29)
