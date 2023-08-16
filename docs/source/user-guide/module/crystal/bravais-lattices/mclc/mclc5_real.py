import radtools as rad

l = rad.lattice_example("MCLC5")
l.plot("primitive", label="primitive")
l.legend()
l.plot("conventional", label="conventional", colour="black")
l.legend()
# Save an image:
l.savefig("mclc5_real.png", elev=20, azim=56, dpi=300, bbox_inches="tight")
# Interactive plot:
l.show(elev=20, azim=56)
