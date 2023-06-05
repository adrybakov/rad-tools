import radtools as rad

l = rad.lattice_example(f"MCLC2")
l.plot("brillouin-kpath")
# Save an image:
l.savefig(
    "mclc2_brillouin.png",
    elev=11,
    azim=64,
    dpi=300,
   bbox_inches="tight",
)
# Interactive plot:
l.show(elev=11, azim=64)
