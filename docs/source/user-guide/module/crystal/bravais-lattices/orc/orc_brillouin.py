import radtools as rad

l = rad.lattice_example(f"ORC")
l.plot("brillouin-kpath")
# Save an image:
l.savefig(
    "orc_brillouin.png",
    elev=35,
    azim=34,
    dpi=300,
   bbox_inches="tight",
)
# Interactive plot:
l.show(elev=35, azim=34)
