import radtools as rad

l = rad.lattice_example(f"ORCF2")
l.plot("brillouin-kpath")
# Save an image:
l.savefig(
    "orcf2_brillouin.png",
    elev=15,
    azim=36,
    dpi=300,
   bbox_inches="tight",
)
# Interactive plot:
l.show(elev=15, azim=36)
