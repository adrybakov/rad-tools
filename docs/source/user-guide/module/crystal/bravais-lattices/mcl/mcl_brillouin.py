import radtools as rad

l = rad.lattice_example(f"MCL")
l.plot("brillouin-kpath")
# Save an image:
l.savefig(
    "mcl_brillouin.png",
    elev=12,
    azim=25,
    dpi=300,
   bbox_inches="tight",
)
# Interactive plot:
l.show(elev=12, azim=25)
