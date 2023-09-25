import radtools as rad

l = rad.lattice_example("ORCF3")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="wigner-seitz")
# Save an image:
backend.save("orcf3_wigner-seitz.png", elev=23, azim=38, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=23, azim=38)
