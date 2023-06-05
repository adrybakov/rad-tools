import radtools as rad

l = rad.lattice_example(f"ORCF3")
l.plot("brillouin-kpath")
# Save an image:
l.savefig(
    "orcf3_brillouin.png",
    elev=25,
    azim=62,
    dpi=300,
   bbox_inches="tight",
)
# Interactive plot:
l.show(elev=25, azim=62)
