import radtools as rad

l = rad.lattice_example("BCC")
backend = rad.PlotlyBackend()
backend.plot(l, kind="brillouin-kpath")
# Save an image:
backend.save("bcc_brillouin.png")
# Interactive plot:
backend.show()
