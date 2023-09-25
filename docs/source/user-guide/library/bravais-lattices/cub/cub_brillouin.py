import radtools as rad

l = rad.lattice_example("CUB")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="brillouin-kpath")
# Save an image:
backend.save("cub_brillouin.png", elev=28, azim=23, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=28, azim=23)
