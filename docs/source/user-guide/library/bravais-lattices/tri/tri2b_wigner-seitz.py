import radtools as rad

l = rad.lattice_example("TRI2b")
backend = rad.PlotlyBackend()
backend.plot(l, kind="wigner-seitz")
# Save an image:
backend.save("tri2b_wigner-seitz.png")
# Interactive plot:
backend.show()
