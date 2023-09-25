import radtools as rad

l = rad.lattice_example("HEX")
backend = rad.PlotlyBackend()
backend.plot(l, kind="primitive")
# Save an image:
backend.save("hex_real.png")
# Interactive plot:
backend.show()
