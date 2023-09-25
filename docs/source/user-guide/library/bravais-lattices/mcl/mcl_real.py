import radtools as rad

l = rad.lattice_example("MCL")
backend = rad.PlotlyBackend()
backend.plot(l, kind="primitive")
# Save an image:
backend.save("mcl_real.png")
# Interactive plot:
backend.show()
