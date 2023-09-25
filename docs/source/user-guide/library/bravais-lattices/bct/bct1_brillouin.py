import radtools as rad

l = rad.lattice_example("BCT1")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="brillouin-kpath")
# Save an image:
backend.save("bct1_brillouin.png", elev=30, azim=28, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=30, azim=28)
