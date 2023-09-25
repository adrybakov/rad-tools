import radtools as rad

l = rad.lattice_example("TRI1a")
backend = rad.PlotlyBackend()
backend.plot(l, kind="primitive")
# Save an image:
backend.save("tri1a_real.png")
# Interactive plot:
backend.show()
