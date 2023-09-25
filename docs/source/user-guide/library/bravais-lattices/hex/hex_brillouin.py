import radtools as rad

l = rad.lattice_example("HEX")
backend = rad.PlotlyBackend()
backend.plot(l, kind="brillouin-kpath")
# Save an image:
backend.save("hex_brillouin.png")
# Interactive plot:
backend.show()
