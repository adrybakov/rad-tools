import radtools as rad

l = rad.lattice_example("RHL2")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="brillouin-kpath")
# Save an image:
backend.save("rhl2_brillouin.png", elev=14, azim=-85, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=14, azim=-85)
