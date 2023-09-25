import radtools as rad

l = rad.lattice_example("ORC")
backend = rad.PlotlyBackend()
backend.plot(l, kind="primitive")
# Save an image:
backend.save("orc_real.png")
# Interactive plot:
backend.show()
