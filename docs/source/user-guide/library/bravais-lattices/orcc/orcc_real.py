import radtools as rad

l = rad.lattice_example("ORCC")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="primitive", label="primitive")
backend.legend()
backend.plot(l, kind="conventional", label="conventional", color="black")
backend.legend()
# Save an image:
backend.save("orcc_real.png", elev=39, azim=15, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=39, azim=15)
