import radtools as rad

l = rad.lattice_example("TRI2b")
backend = rad.PlotlyBackend()
backend.plot(l, kind="brillouin-kpath")
# Save an image:
backend.save("tri2b_brillouin.png")
# Interactive plot:
backend.show()
