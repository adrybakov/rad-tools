import radtools as rad

l = rad.lattice_example("RHL2")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="primitive")
# Save an image:
backend.save("rhl2_real.png", elev=35, azim=52, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=35, azim=52)
