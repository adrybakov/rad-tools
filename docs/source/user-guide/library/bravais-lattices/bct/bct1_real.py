import radtools as rad

l = rad.lattice_example("BCT1")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="primitive", label="primitive")
backend.legend()
backend.plot(l, kind="conventional", label="conventional", color="black")
backend.legend()
# Save an image:
backend.save("bct1_real.png", elev=37, azim=72, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=37, azim=72)
