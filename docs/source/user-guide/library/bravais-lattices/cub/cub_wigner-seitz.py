import radtools as rad

l = rad.lattice_example("CUB")
backend = rad.PlotlyBackend()
backend.plot(l, kind="wigner-seitz")
# Save an image:
backend.save("cub_wigner-seitz.png")
# Interactive plot:
backend.show()
