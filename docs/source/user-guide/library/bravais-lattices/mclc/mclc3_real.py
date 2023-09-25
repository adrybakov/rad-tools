import radtools as rad

l = rad.lattice_example("MCLC3")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="primitive", label="primitive")
backend.legend()
backend.plot(l, kind="conventional", label="conventional", color="black")
backend.legend()
# Save an image:
backend.save("mclc3_real.png", elev=39, azim=68, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=39, azim=68)
