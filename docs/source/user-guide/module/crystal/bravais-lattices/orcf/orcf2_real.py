import radtools as rad

l = rad.lattice_example("ORCF2")
l.plot("primitive", label="primitive")
l.legend()
l.plot("conventional", label="conventional", colour="black")
l.legend()
# Save an image:
l.savefig("orcf2_real.png", elev=25, azim=28, dpi=300, bbox_inches="tight")
# Interactive plot:
l.show(elev=25, azim=28)
