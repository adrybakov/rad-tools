import radtools as rad

l = rad.lattice_example("MCLC5")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="primitive", label="primitive")
backend.legend()
backend.plot(l, kind="conventional", label="conventional", color="black")
backend.legend()
# Save an image:
backend.save("mclc5_real.png", elev=20, azim=56, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=20, azim=56)
