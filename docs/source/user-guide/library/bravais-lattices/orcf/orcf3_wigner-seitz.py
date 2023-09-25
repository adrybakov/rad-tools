import radtools as rad

l = rad.lattice_example("ORCF3")
backend = rad.PlotlyBackend()
backend.plot(l, kind="wigner-seitz")
# Save an image:
backend.save("orcf3_wigner-seitz.png")
# Interactive plot:
backend.show()
