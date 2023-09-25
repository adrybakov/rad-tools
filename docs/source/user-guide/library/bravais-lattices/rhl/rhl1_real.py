import radtools as rad

l = rad.lattice_example("RHL1")
backend = rad.PlotlyBackend()
backend.plot(l, kind="primitive")
# Save an image:
backend.save("rhl1_real.png")
# Interactive plot:
backend.show()
