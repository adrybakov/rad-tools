import radtools as rad

l = rad.lattice_example("ORCF2")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="wigner-seitz")
# Save an image:
backend.save("orcf2_wigner-seitz.png", elev=38, azim=14, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=38, azim=14)
