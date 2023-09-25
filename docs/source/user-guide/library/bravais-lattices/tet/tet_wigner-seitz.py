import radtools as rad

l = rad.lattice_example("TET")
backend = rad.PlotlyBackend()
backend.plot(l, kind="wigner-seitz")
# Save an image:
backend.save("tet_wigner-seitz.png")
# Interactive plot:
backend.show()
