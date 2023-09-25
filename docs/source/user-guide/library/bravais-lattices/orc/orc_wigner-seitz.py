import radtools as rad

l = rad.lattice_example("ORC")
backend = rad.PlotlyBackend()
backend.plot(l, kind="wigner-seitz")
# Save an image:
backend.save("orc_wigner-seitz.png")
# Interactive plot:
backend.show()
