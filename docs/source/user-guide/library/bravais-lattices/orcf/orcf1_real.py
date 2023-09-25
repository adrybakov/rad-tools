import radtools as rad

l = rad.lattice_example("ORCF1")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="primitive", label="primitive")
backend.legend()
backend.plot(l, kind="conventional", label="conventional", color="black")
backend.legend()
# Save an image:
backend.save("orcf1_real.png", elev=24, azim=38, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=24, azim=38)
