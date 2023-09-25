import radtools as rad

l = rad.lattice_example("ORC")
backend = rad.PlotlyBackend()
backend.plot(l, kind="brillouin-kpath")
# Save an image:
backend.save("orc_brillouin.png")
# Interactive plot:
backend.show()
