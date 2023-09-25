import radtools as rad

l = rad.lattice_example("CUB")
backend = rad.PlotlyBackend()
backend.plot(l, kind="primitive")
# Save an image:
backend.save("cub_real.png")
# Interactive plot:
backend.show()
