import radtools as rad

l = rad.lattice_example("MCLC4")
backend = rad.PlotlyBackend()
backend.plot(l, kind="brillouin-kpath")
# Save an image:
backend.save("mclc4_brillouin.png")
# Interactive plot:
backend.show()
