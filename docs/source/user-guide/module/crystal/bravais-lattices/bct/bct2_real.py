import radtools as rad

l = rad.lattice_example("BCT2")
l.plot("primitive", label="primitive")
l.legend()
l.plot("conventional", label="conventional", colour="black")
l.legend()
# Save an image:
l.savefig("bct2_real.png", elev=40, azim=85, dpi=300, bbox_inches="tight")
# Interactive plot:
l.show(elev=40, azim=85)
