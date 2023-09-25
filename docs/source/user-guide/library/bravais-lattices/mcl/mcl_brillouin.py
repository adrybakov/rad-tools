import radtools as rad

l = rad.lattice_example("MCL")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="brillouin-kpath")
# Save an image:
backend.save("mcl_brillouin.png", elev=12, azim=25, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=12, azim=25)
