import radtools as rad

l = rad.lattice_example("ORCF2")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="primitive", label="primitive")
backend.legend()
backend.plot(l, kind="conventional", label="conventional", color="black")
backend.legend()
# Save an image:
backend.save("orcf2_real.png", elev=25, azim=28, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=25, azim=28)
