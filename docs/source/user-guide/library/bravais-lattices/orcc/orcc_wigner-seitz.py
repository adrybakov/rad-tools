import radtools as rad

l = rad.lattice_example("ORCC")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="wigner-seitz")
# Save an image:
backend.save("orcc_wigner-seitz.png", elev=33, azim=-40, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=33, azim=-40)
