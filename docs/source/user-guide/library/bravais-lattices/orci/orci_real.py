import radtools as rad

l = rad.lattice_example("ORCI")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="primitive", label="primitive")
backend.legend()
backend.plot(l, kind="conventional", label="conventional", color="black")
backend.legend()
# Save an image:
backend.save("orci_real.png", elev=32, azim=-12, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=32, azim=-12)
