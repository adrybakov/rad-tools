import radtools as rad

l = rad.lattice_example("FCC")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="primitive", label="primitive")
backend.legend()
backend.plot(l, kind="conventional", label="conventional", color="black")
backend.legend()
# Save an image:
backend.save("fcc_real.png", elev=28, azim=23, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=28, azim=23)
