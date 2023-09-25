import radtools as rad

l = rad.lattice_example("TET")
backend = rad.PlotlyBackend()
backend.plot(l, kind="brillouin-kpath")
# Save an image:
backend.save("tet_brillouin.png")
# Interactive plot:
backend.show()
