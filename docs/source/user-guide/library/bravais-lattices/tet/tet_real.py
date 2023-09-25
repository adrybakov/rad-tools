import radtools as rad

l = rad.lattice_example("TET")
backend = rad.PlotlyBackend()
backend.plot(l, kind="primitive")
# Save an image:
backend.save("tet_real.png")
# Interactive plot:
backend.show()
