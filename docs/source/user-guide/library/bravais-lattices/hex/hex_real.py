import radtools as rad

l = rad.lattice_example("HEX")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="primitive")
# Save an image:
backend.save("hex_real.png", elev=35, azim=23, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=35, azim=23)
