import radtools as rad

l = rad.lattice_example("MCL")
backend = rad.PlotlyBackend()
backend.plot(l, kind="wigner-seitz")
# Save an image:
backend.save("mcl_wigner-seitz.png")
# Interactive plot:
backend.show()
