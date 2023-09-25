import radtools as rad

l = rad.lattice_example("MCL")
backend = rad.PlotlyBackend()
backend.plot(l, kind="brillouin-kpath")
# Save an image:
backend.save("mcl_brillouin.png")
# Interactive plot:
backend.show()
