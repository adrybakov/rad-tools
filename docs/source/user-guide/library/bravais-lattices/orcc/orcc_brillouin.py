import radtools as rad

l = rad.lattice_example("ORCC")
backend = rad.PlotlyBackend()
backend.plot(l, kind="brillouin-kpath")
# Save an image:
backend.save("orcc_brillouin.png")
# Interactive plot:
backend.show()
