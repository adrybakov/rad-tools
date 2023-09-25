import radtools as rad

l = rad.lattice_example("BCT2")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="brillouin-kpath")
# Save an image:
backend.save("bct2_brillouin.png", elev=36, azim=28, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=36, azim=28)
