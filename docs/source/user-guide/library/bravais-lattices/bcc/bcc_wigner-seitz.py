import radtools as rad

l = rad.lattice_example("BCC")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="wigner-seitz")
# Save an image:
backend.save("bcc_wigner-seitz.png", elev=46, azim=19, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=46, azim=19)
