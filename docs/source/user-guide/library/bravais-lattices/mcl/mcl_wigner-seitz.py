import radtools as rad

l = rad.lattice_example("MCL")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="wigner-seitz")
# Save an image:
backend.save("mcl_wigner-seitz.png", elev=11, azim=37, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=11, azim=37)
