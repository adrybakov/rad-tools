import radtools as rad

l = rad.lattice_example("MCLC5")
backend = rad.PlotlyBackend()
backend.plot(l, kind="wigner-seitz")
# Save an image:
backend.save("mclc5_wigner-seitz.png")
# Interactive plot:
backend.show()
