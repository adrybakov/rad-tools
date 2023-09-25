import radtools as rad

l = rad.lattice_example("MCLC2")
backend = rad.PlotlyBackend()
backend.plot(l, kind="brillouin-kpath")
# Save an image:
backend.save("mclc2_brillouin.png")
# Interactive plot:
backend.show()
