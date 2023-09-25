import radtools as rad

l = rad.lattice_example("TRI1b")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="primitive")
# Save an image:
backend.save("tri1b_real.png", elev=12, azim=11, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=12, azim=11)
