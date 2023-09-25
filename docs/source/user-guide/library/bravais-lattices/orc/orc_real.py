import radtools as rad

l = rad.lattice_example("ORC")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="primitive")
# Save an image:
backend.save("orc_real.png", elev=36, azim=35, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=36, azim=35)
