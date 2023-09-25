import radtools as rad

l = rad.lattice_example("FCC")
backend = rad.PlotlyBackend()
backend.plot(l, kind="wigner-seitz")
# Save an image:
backend.save("fcc_wigner-seitz.png")
# Interactive plot:
backend.show()
