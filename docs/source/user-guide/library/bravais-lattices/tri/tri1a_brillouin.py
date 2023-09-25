import radtools as rad

l = rad.lattice_example("TRI1a")
backend = rad.PlotlyBackend()
backend.plot(l, kind="brillouin-kpath")
# Save an image:
backend.save("tri1a_brillouin.png")
# Interactive plot:
backend.show()
