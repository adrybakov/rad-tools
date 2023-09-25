import radtools as rad

l = rad.lattice_example("TRI1a")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="brillouin-kpath")
# Save an image:
backend.save("tri1a_brillouin.png", elev=45, azim=-55, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=45, azim=-55)
