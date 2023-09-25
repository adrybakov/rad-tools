import radtools as rad

l = rad.lattice_example("TRI2a")
backend = rad.PlotlyBackend()
backend.plot(l, kind="brillouin-kpath")
# Save an image:
backend.save("tri2a_brillouin.png")
# Interactive plot:
backend.show()
