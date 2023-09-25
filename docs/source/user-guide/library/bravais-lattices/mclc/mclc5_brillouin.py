import radtools as rad

l = rad.lattice_example("MCLC5")
backend = rad.PlotlyBackend()
backend.plot(l, kind="brillouin-kpath")
# Save an image:
backend.save("mclc5_brillouin.png")
# Interactive plot:
backend.show()
