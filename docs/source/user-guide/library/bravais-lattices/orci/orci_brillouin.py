import radtools as rad

l = rad.lattice_example("ORCI")
backend = rad.PlotlyBackend()
backend.plot(l, kind="brillouin-kpath")
# Save an image:
backend.save("orci_brillouin.png")
# Interactive plot:
backend.show()
