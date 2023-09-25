import radtools as rad

l = rad.lattice_example("TRI1b")
backend = rad.PlotlyBackend()
backend.plot(l, kind="primitive")
# Save an image:
backend.save("tri1b_real.png")
# Interactive plot:
backend.show()
