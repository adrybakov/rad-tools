import radtools as rad

l = rad.lattice_example("CUB")
backend = rad.PlotlyBackend()
backend.plot(l, kind="brillouin-kpath")
# Save an image:
backend.save("cub_brillouin.png")
# Interactive plot:
backend.show()
