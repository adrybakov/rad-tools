import radtools as rad

l = rad.lattice_example("TRI2a")
backend = rad.PlotlyBackend()
backend.plot(l, kind="wigner-seitz")
# Save an image:
backend.save("tri2a_wigner-seitz.png")
# Interactive plot:
backend.show()
