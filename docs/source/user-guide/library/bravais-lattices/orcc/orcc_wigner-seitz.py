import radtools as rad

l = rad.lattice_example("ORCC")
backend = rad.PlotlyBackend()
backend.plot(l, kind="wigner-seitz")
# Save an image:
backend.save("orcc_wigner-seitz.png")
# Interactive plot:
backend.show()
