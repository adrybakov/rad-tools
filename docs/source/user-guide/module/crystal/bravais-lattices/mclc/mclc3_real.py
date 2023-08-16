import radtools as rad

l = rad.lattice_example("MCLC3")
l.plot("primitive", label="primitive")
l.legend()
l.plot("conventional", label="conventional", colour="black")
l.legend()
# Save an image:
l.savefig("mclc3_real.png", elev=39, azim=68, dpi=300, bbox_inches="tight")
# Interactive plot:
l.show(elev=39, azim=68)
