import radtools as rad

l = rad.lattice_example("FCC")
backend = rad.PlotlyBackend()
backend.plot(l, kind="brillouin-kpath")
# Save an image:
backend.save("fcc_brillouin.png")
# Interactive plot:
backend.show()
