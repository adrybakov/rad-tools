import radtools as rad

l = rad.lattice_example("ORCF1")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="wigner-seitz")
# Save an image:
backend.save("orcf1_wigner-seitz.png", elev=44, azim=28, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=44, azim=28)
