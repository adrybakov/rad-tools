import radtools as rad

l = rad.lattice_example("MCLC1")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="primitive", label="primitive")
backend.legend()
backend.plot(l, kind="conventional", label="conventional", color="black")
backend.legend()
# Save an image:
backend.save("mclc1_real.png", elev=29, azim=83, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=29, azim=83)
