import radtools as rad

l = rad.lattice_example("MCLC4")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="brillouin-kpath")
# Save an image:
backend.save("mclc4_brillouin.png", elev=19, azim=45, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=19, azim=45)
