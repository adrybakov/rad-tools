import radtools as rad

l = rad.lattice_example("ORCI")
backend = rad.PlotlyBackend()
backend.plot(l, kind="wigner-seitz")
# Save an image:
backend.save("orci_wigner-seitz.png")
# Interactive plot:
backend.show()
