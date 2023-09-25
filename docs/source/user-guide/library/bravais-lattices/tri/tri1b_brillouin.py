import radtools as rad

l = rad.lattice_example("TRI1b")
backend = rad.PlotlyBackend()
backend.plot(l, kind="brillouin-kpath")
# Save an image:
backend.save("tri1b_brillouin.png")
# Interactive plot:
backend.show()
