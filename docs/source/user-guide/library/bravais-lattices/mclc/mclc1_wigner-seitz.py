import radtools as rad

l = rad.lattice_example("MCLC1")
backend = rad.PlotlyBackend()
backend.plot(l, kind="wigner-seitz")
# Save an image:
backend.save("mclc1_wigner-seitz.png")
# Interactive plot:
backend.show()
