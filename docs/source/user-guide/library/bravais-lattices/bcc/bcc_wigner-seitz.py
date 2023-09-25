import radtools as rad

l = rad.lattice_example("BCC")
backend = rad.PlotlyBackend()
backend.plot(l, kind="wigner-seitz")
# Save an image:
backend.save("bcc_wigner-seitz.png")
# Interactive plot:
backend.show()
