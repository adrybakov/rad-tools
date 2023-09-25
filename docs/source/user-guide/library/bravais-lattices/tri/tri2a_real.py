import radtools as rad

l = rad.lattice_example("TRI2a")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="primitive")
# Save an image:
backend.save("tri2a_real.png", elev=39, azim=44, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=39, azim=44)
