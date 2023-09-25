import radtools as rad

l = rad.lattice_example("MCLC3")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="wigner-seitz")
# Save an image:
backend.save("mclc3_wigner-seitz.png", elev=41, azim=61, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=41, azim=61)
