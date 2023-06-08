import radtools as rad

l = rad.lattice_example(f"BCT1")
l.plot("brillouin-kpath")
# Save an image:
l.savefig(
    "bct1_brillouin.png",
    elev=30,
    azim=28,
    dpi=300,
   bbox_inches="tight",
)
# Interactive plot:
l.show(elev=30, azim=28)
