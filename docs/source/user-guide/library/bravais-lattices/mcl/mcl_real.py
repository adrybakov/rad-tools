import radtools as rad

l = rad.lattice_example("MCL")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="primitive")
# Save an image:
backend.save("mcl_real.png", elev=25, azim=40, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=25, azim=40)
