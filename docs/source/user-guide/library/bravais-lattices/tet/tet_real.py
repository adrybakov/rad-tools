import radtools as rad

l = rad.lattice_example("TET")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="primitive")
# Save an image:
backend.save("tet_real.png", elev=30, azim=30, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=30, azim=30)
