import radtools as rad

l = rad.lattice_example("RHL1")
backend = rad.PlotlyBackend()
backend.plot(l, kind="brillouin-kpath")
# Save an image:
backend.save("rhl1_brillouin.png")
# Interactive plot:
backend.show()
