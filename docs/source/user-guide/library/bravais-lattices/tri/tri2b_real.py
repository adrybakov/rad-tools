import radtools as rad

l = rad.lattice_example("TRI2b")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="primitive")
# Save an image:
backend.save("tri2b_real.png", elev=17, azim=54, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=17, azim=54)
