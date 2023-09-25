import radtools as rad

l = rad.lattice_example("TRI2b")
backend = rad.PlotlyBackend()
backend.plot(l, kind="primitive")
# Save an image:
backend.save("tri2b_real.png")
# Interactive plot:
backend.show()
