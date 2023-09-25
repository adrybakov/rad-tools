import radtools as rad

l = rad.lattice_example("RHL2")
backend = rad.PlotlyBackend()
backend.plot(l, kind="primitive")
# Save an image:
backend.save("rhl2_real.png")
# Interactive plot:
backend.show()
