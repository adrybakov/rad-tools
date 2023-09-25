import radtools as rad

l = rad.lattice_example("RHL1")
backend = rad.PlotlyBackend()
backend.plot(l, kind="wigner-seitz")
# Save an image:
backend.save("rhl1_wigner-seitz.png")
# Interactive plot:
backend.show()
