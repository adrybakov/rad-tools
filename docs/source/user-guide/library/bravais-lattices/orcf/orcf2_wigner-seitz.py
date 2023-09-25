import radtools as rad

l = rad.lattice_example("ORCF2")
backend = rad.PlotlyBackend()
backend.plot(l, kind="wigner-seitz")
# Save an image:
backend.save("orcf2_wigner-seitz.png")
# Interactive plot:
backend.show()
