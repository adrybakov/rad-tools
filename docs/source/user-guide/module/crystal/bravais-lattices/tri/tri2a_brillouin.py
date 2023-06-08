import radtools as rad

l = rad.lattice_example(f"TRI2a")
l.plot("brillouin-kpath")
# Save an image:
l.savefig(
    "tri2a_brillouin.png",
    elev=10,
    azim=32,
    dpi=300,
   bbox_inches="tight",
)
# Interactive plot:
l.show(elev=10, azim=32)
