import radtools as rad

l = rad.lattice_example("TRI1b")
backend = rad.PlotlyBackend()
backend.plot(l, kind="wigner-seitz")
# Save an image:
backend.save("tri1b_wigner-seitz.png")
# Interactive plot:
backend.show()
