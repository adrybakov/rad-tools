import radtools as rad

l = rad.lattice_example("MCLC4")
backend = rad.PlotlyBackend()
backend.plot(l, kind="wigner-seitz")
# Save an image:
backend.save("mclc4_wigner-seitz.png")
# Interactive plot:
backend.show()
