import radtools as rad

l = rad.lattice_example(f"FCC")
l.plot("brillouin-kpath")
# Save an image:
l.savefig(
    "fcc_brillouin.png",
    elev=23,
    azim=28,
    dpi=300,
   bbox_inches="tight",
)
# Interactive plot:
l.show(elev=23, azim=28)
