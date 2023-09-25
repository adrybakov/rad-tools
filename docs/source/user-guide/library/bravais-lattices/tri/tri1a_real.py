import radtools as rad

l = rad.lattice_example("TRI1a")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="primitive")
# Save an image:
backend.save("tri1a_real.png", elev=31, azim=-20, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=31, azim=-20)
