import radtools as rad

l = rad.lattice_example("ORCF1")
l.plot("primitive", label="primitive")
l.legend()
l.plot("conventional", label="conventional", colour="black")
l.legend()
# Save an image:
l.savefig("orcf1_real.png", elev=24, azim=38, dpi=300, bbox_inches="tight")
# Interactive plot:
l.show(elev=24, azim=38)
