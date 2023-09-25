import radtools as rad

l = rad.lattice_example("ORCF3")
backend = rad.PlotlyBackend()
backend.plot(l, kind="brillouin-kpath")
# Save an image:
backend.save("orcf3_brillouin.png")
# Interactive plot:
backend.show()
