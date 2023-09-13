import radtools as rad

l = rad.lattice_example("ORCI")
l.plot("primitive", label="primitive")
l.legend()
l.plot("conventional", label="conventional", colour="black")
l.legend()
# Save an image:
l.savefig("orci_real.png", elev=32, azim=-12, dpi=300, bbox_inches="tight")
# Interactive plot:
l.show(elev=32, azim=-12)
