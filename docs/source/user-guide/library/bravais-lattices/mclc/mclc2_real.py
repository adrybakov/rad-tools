import radtools as rad

l = rad.lattice_example("MCLC2")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="primitive", label="primitive")
backend.legend()
backend.plot(l, kind="conventional", label="conventional", color="black")
backend.legend()
# Save an image:
backend.save("mclc2_real.png", elev=30, azim=62, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=30, azim=62)
