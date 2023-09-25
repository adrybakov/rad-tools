import radtools as rad

l = rad.lattice_example("MCLC1")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="wigner-seitz")
# Save an image:
backend.save("mclc1_wigner-seitz.png", elev=45, azim=67, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=45, azim=67)
