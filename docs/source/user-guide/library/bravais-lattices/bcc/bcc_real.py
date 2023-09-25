import radtools as rad

l = rad.lattice_example("BCC")
backend = rad.MatplotlibBackend()
backend.plot(l, kind="primitive", label="primitive")
backend.legend()
backend.plot(l, kind="conventional", label="conventional", color="black")
backend.legend()
# Save an image:
backend.save("bcc_real.png", elev=30, azim=-180, dpi=300, bbox_inches="tight")
# Interactive plot:
backend.show(elev=30, azim=-180)
